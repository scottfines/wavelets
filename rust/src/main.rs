mod haar;
mod arrays;

fn main() {
    let data = [1,3,5,11,12,13,0,1];

    let wavelet = haar::cascade(&data);
    wavelet.iter().for_each(|k| { print!("{},",k) });
    println!("");

    let inverse = haar::inverse_cascade(&wavelet);
    inverse.iter().for_each(|k| { print!("{},",k) });
    println!("");
}
