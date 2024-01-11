# Image quilting for texture synthesis and transfer

This repository is an implementation of the image quilting for texture synthesis and transfer using Cpp.

## Reference:

1. Alexei A. Efros and William T. Freeman. 2001. Image quilting for texture synthesis and transfer. In Proceedings of the 28th annual conference on Computer graphics and interactive techniques (SIGGRAPH '01). Association for Computing Machinery, New York, NY, USA, 341–346. https://doi.org/10.1145/383259.383296

2. https://github.com/rohitrango/Image-Quilting-for-Texture-Synthesis

## Development Environment

MacOS 10.15

## How to run:

```
make all
make run [ARGS="[-i input_path] [-o output_path] [-p patch_size] [-v overlap_ratio] [-s scale] [-t tolerance] [-n n_iter]"]
```

`[ ]` means optional, you can run through `make run`, `make run ARGS="-p 20"`, `make run ARGS="-p 20 -s 5"`. Anyone of them is runable.

## Argument for Executable:

- `-i` for input image path.
- `-o` for output image path.

## Simple Example:

Input image:

<img class="myImage" src="input/case1.png" /> <img class="myImage" src="input/case2.png" /> <img class="myImage" src="input/case3.png" />

Output synthesized image:

<img class="myImage" src="output/case1.png" /> <img class="myImage" src="output/case2.png" /> <img class="myImage" src="output/case3.png" />

## Note

The effect may varies with respect to our definition on patch size, error tolrence and max iteration. I haven't tested proper parameters for each case. Maybe you can try on different parameters for better results.
