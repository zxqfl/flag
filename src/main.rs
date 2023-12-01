use rayon::prelude::*;
use std::sync::Mutex;
use std::f64::consts::PI;
use image::{RgbImage, Rgb};
use std::collections::HashMap;

const VIDEO: bool = false;

struct PseudoRandomRounder {
    n: u64,
}

impl PseudoRandomRounder {
    fn new() -> Self {
        PseudoRandomRounder {
            n: 0,
        }
    }

    fn next(&mut self) -> f64 {
        let mut result = 0.0;
        let mut x = self.n;
        let mut amt = 1.;
        self.n += 1;
        while x > 0 {
            amt /= 2.;
            if x % 2 != 0 {
                result += amt;
            }
            x /= 2;
        }
        result
    }

    fn round(&mut self, x: f64) -> u32 {
        let x = x + self.next();
        let x = x.floor();
        assert!(x >= 0.);
        x as u32
    }
}

fn rem(a: i64, b: i64) -> i64 {
    let mut x = a % b;
    if x < 0 {
        x += b;
    }
    x
}

fn render_star(sz: u32) -> tiny_skia::Pixmap {
    let opt = usvg::Options::default();
    let svg_data: &[u8] = include_bytes!("./star.svg");
    let rtree = usvg::Tree::from_data(&svg_data, &opt.to_ref()).unwrap();
    let mut pixmap = tiny_skia::Pixmap::new(sz, sz).unwrap();
    resvg::render(&rtree, usvg::FitTo::Size(sz, sz), pixmap.as_mut()).unwrap();
    pixmap
}

fn fmin(a: f64, b: f64) -> f64 {
    if a < b {
        a
    } else {
        b
    }
}
fn fmax(a: f64, b: f64) -> f64 {
    -fmin(-a, -b)
}

fn dist2(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    let dx = x1 - x2;
    let dy = y1 - y2;
    dx * dx + dy * dy
}

fn dist(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    dist2(x1, y1, x2, y2).sqrt()
}

#[derive(Clone, Copy, Debug)]
struct Star {
    x: f64,
    y: f64,
    radius: f64,
}

type Point = (f64, f64);

struct CircleWalker {
    centre1: Point,
    radius1: f64,
    centre2: Point,
    radius2: f64,
    angle1: f64,
    angle2: f64,
}

impl CircleWalker {
    fn new(centre1: Point, radius1: f64, centre2: Point, radius2: f64) -> Self {
        CircleWalker {
            centre1, centre2, radius1, radius2,
            angle1: 0.0,
            angle2: 0.0,
        }
    }

    fn set_angle(&mut self, angle1: f64) {
        self.angle1 = angle1;
        // self.angle2 = angle1;
    }

    fn get_points(&mut self) -> (Point, Point) {
        let (cx1, cy1) = self.centre1;
        let (cx2, cy2) = self.centre2;
        let x1 = cx1 + self.angle1.cos() * self.radius1;
        let y1 = cy1 + self.angle1.sin() * self.radius1;
        let mut step = 1.;
        while step > 1e-10 {
            let x2 = cx2 + self.angle2.cos() * self.radius2;
            let y2 = cy2 + self.angle2.sin() * self.radius2;
            let dx = -self.angle2.sin() * self.radius2;
            let dy = self.angle2.cos() * self.radius2;
            let inc_or_dec = dx * (x1 - x2) + dy * (y1 - y2) > 0.;
            let new_angle2;
            if inc_or_dec {
                new_angle2 = self.angle2 + step;
            } else {
                new_angle2 = self.angle2 - step;
            }
            let new_x2 = cx2 + new_angle2.cos() * self.radius2;
            let new_y2 = cy2 + new_angle2.sin() * self.radius2;
            let d_new = dist2(x1, y1, new_x2, new_y2);
            let d_old = dist2(x1, y1, x2, y2);
            if d_new < d_old {
                self.angle2 = new_angle2;
            } else {
                step *= 0.3;
            }
        }
        let x2 = cx2 + self.angle2.cos() * self.radius2;
        let y2 = cy2 + self.angle2.sin() * self.radius2;
        ((x1, y1), (x2, y2))
    }
}

fn blend(mut a: Rgb<u8>, b: Rgb<u8>, alpha: u8) -> Rgb<u8> {
    for i in 0..3 {
        let diff = (b.0[i] as i32 - a.0[i] as i32) * alpha as i32 / 255;
        a.0[i] = (a.0[i] as i32 + diff) as u8;
    }
    a
}

fn main() {
    let integral_steps = 40;
    let downsample = 4;
    let scale = 32 * downsample;
    let w = 13 * 19 * 2 * scale;
    let h = 13 * 10 * 2 * scale;
    let union_h = h * 7 / 13;
    let union_w = w * 4 / 10;

    let blue = Rgb([0, 0x28, 0x68]);
    // let blue_slightly_white = Rgb([0x05, 0x2c, 0x6c]);
    // let blue_slightly_white = Rgb([0x10, 0x48, 0x88]);
    let blue_slightly_white = blue;
    let blue_whitened = Rgb([0x15, 0x38, 0x78]);

    let red = Rgb([0xbf, 0x0a, 0x30]);
    let white = Rgb([255, 255, 255]);

    let mut filenames = Vec::new();

    let mut img = RgbImage::new(w, h);

    const KCIRC: f64 = 3000.;
    for x in 0..w {
        if x % 100 == 0 {
            println!("{}/{}", x+1, w);
        }
        for y in 0..h {
            let color;
            let in_union = x < union_w && y < union_h;
            // if in_union {
            //     color = blue;
            // } else {
            {
                let cx = union_w as f64 / 2.;
                let cy = union_h as f64 / 2.;

                let x = x as f64 - cx;
                let y = y as f64 - cy;

                let scale = scale as f64;

                let alpha = (x*x + y*y) / (2.*y);
                // let stripe_width = w as f64 / 100.;
                // let idx = (alpha / stripe_width) as i64;
                let idx = (1. / alpha * scale * KCIRC).round() as i64;
                let n13 = rem(idx + 2, 13);
                if in_union {
                    if n13 == 0 || n13 == 12 {
                        color = blue_slightly_white;
                    } else {
                        color = blue;
                    }
                    // color = blue;
                } else {
                    if n13 % 2 == 0 {
                        color = red;
                    } else {
                        color = white;
                    }
                }
            }
            img.put_pixel(x, y, color);
        }
    }

    if downsample > 1 {
        println!("Downsampling...");
        img = image::imageops::resize(&img, img.width() / downsample, img.height() / downsample, image::imageops::FilterType::Triangle);
    }
    let w = w / downsample;
    let h = h / downsample;
    let union_w = union_w / downsample;
    let union_h = union_h / downsample;
    let scale = scale / downsample;
    let range = if VIDEO {
        0..integral_steps
    } else {
        0..1
    };

    for t in range {
        let mut img = img.clone();

        let scale = scale as f64;
        let nk = 200;
        let done = Mutex::new(0);
        let stars: Vec<Vec<Star>> = (-nk..=nk).into_par_iter()
            .map(|k| {
                let mut stars: Vec<Star> = Vec::new();
                {
                    let mut done = done.lock().unwrap();
                    *done += 1;
                    println!("{done}/{}", 2*nk+1);
                }
                let n1 = (k*13) as f64 + 12.;
                let n2 = (k*13) as f64 + 13.;
                let r1 = 1. / (n1 - 2.5) * scale * KCIRC;
                let r2 = 1. / (n2 - 1.5) * scale * KCIRC;
                let cx = union_w as f64 / 2.;
                let cy = union_h as f64 / 2.;
                for dir in [-1., 1.] {
                    let mut walker = CircleWalker::new(
                        (cx, cy + r1),
                        r1,
                        (cx, cy + r2),
                        r2,
                    );
                    let start_angle = -PI / 2.;
                    let smallstep_mul = 100.;
                    // let mut angle = start_angle + 0.1 * dir;
                    let mut angle = start_angle;
                    for it in 0.. {
                        walker.set_angle(angle);
                        let ((x1, y1), (x2, y2)) = walker.get_points();
                        let x = (x1 + x2) / 2.;
                        let y = (y1 + y2) / 2.;
                        let radius = dist(x1, y1, x2, y2);
                        if (it+t) % integral_steps == integral_steps / 2 || it == 0 {
                            if radius < union_w as f64 {
                                let radius = radius * 1.5;
                                if x+radius >= 0. && y+radius >= 0. && x-radius < union_w as f64 && y-radius < union_h as f64 {
                                    if (x - cx).abs() > radius * 1.4 {
                                        stars.push(Star { x, y, radius });
                                    }
                                    // else {
                                    //     stars.push(Star { x: cx, y, radius });
                                    // }
                                }
                            }
                        }

                        angle += dir * fmax(0.2 / smallstep_mul, radius) / r1.abs() * (8. / integral_steps as f64);

                        if (angle - start_angle) * dir > PI {
                            break;
                        }
                    }
                }
                // let nstar = 320;
                // for i in 0..nstar {
                //     let cx = union_w as f64 / 2.;
                //     let cy = union_h as f64 / 2.;

                //     let frac = i as f64 / nstar as f64;
                //     let angle = frac * PI * 2.;
                //     let x1 = cx + angle.cos() * r1;
                //     let y1 = cy + angle.sin() * r1 + r1;
                //     let x2 = cx + angle.cos() * r2;
                //     let y2 = cy + angle.sin() * r2 + r2;
                //     let x = (x1 + x2) / 2.;
                //     let y = (y1 + y2) / 2.;
                //     if !(x >= 0. && y >= 0. && x < union_w as f64 && y < union_h as f64) {
                //         continue;
                //     }
                //     let radius = 0.1 * fmin(dist(x, y, x1, y1), dist(x, y, x2, y2));
                //     let star = Star { x, y, radius };
                //     stars.push(star);
                //     // println!("r1={:.3} r2={:.3} {:?}", r1, r2, star);
                // }
                stars
            })
            .collect();
        let mut star_cache: HashMap<u32, tiny_skia::Pixmap> = HashMap::new();
        for stars in stars {
            let mut radius_rounder = PseudoRandomRounder::new();
            for Star { x, y, radius } in stars {
                if !VIDEO {
                    let r = radius * 2.;
                    if x < r || y < r || x + r > union_w as f64 || y + r > union_h as f64 {
                        continue;
                    }
                }
                let x = x.round() as i64;
                let y = y.round() as i64;
                // let iradius = radius.round() as i64;

                // let x = x_rounder.round(x);
                // let y = y_rounder.round(y);

                // let iradius = randomized_round(radius);
                // if iradius <= 0 {
                //     continue;
                // }
                // let star_size = (iradius * 2 + 1) as u32;
                let star_size = radius_rounder.round(radius * 2.);
                if star_size == 0 {
                    continue;
                }
                if !star_cache.contains_key(&star_size) {
                    star_cache.insert(star_size, render_star(star_size));
                }
                let star_bitmap = &star_cache[&star_size];
                for star_x in 0..star_size {
                    for star_y in 0..star_size {
                        let wx = x - star_size as i64 / 2 + star_x as i64;
                        let wy = y - star_size as i64 / 2 + star_y as i64;
                        if 0 <= wx && 0 <= wy && wx < union_w as i64 && wy < union_h as i64 {
                            let alpha = star_bitmap.pixel(star_x, star_y).unwrap().alpha();
                            let mut px = *img.get_pixel(wx as u32, wy as u32);
                            px = blend(px, white, alpha);
                            img.put_pixel(wx as u32, wy as u32, px);
                        }
                    }
                }
            }
        }

        for x in 0..union_w {
            for y in 0..union_h {
                let cx = union_w as f64 / 2.;
                let cy = union_h as f64 / 2.;
                let dx = x as f64;
                let dy = y as f64;
                let d = dist(cx, cy, dx, dy);
                let lim = union_w as f64 / 7.;
                if d < lim {
                    let mul = d / lim;
                    let mul = 1. - mul;
                    let mul = mul * mul;
                    let alpha = (mul * 255.) as u8;
                    let mut px = *img.get_pixel(x, y);
                    px = blend(px, blue_whitened, alpha);
                    img.put_pixel(x, y, px);
                }
            }
        }
        let filename;
        if VIDEO {
            filename = format!("test-{:05}.png", t);
        } else {
            filename = "test.png".to_string();
        }
        println!("Saving {}", filename);
        img.save(&filename).unwrap();
        filenames.push(filename);
    }
    for name in filenames {
        println!("{}", name);
    }
}
