# CROC
Crazy Rust Obvious Chemistry

Chemistry Rusty Obvious Crocodile

Cool Rust Open Chemistry 

Basically, it probably means something...


## Developers
Most of the difficulty here arises due to the r-libcint project (a wrapper around the libcint electron integral library). When this is moved into its own crate, development will likely get much easier. 

### Prereqs
cmake: necessary, for now, in order to use libcint for electron integrals.
rust-nightly: this is necessary for [bindgen](https://crates.io/crates/bindgen) (I think)

### Super important git command
In order to get the libcint library you must recusively call the submodules during the clone.
``` bash
git clone --recurse-submodules <this repo>

```


