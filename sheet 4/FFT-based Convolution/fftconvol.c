#include <stdlib.h>
#include <fftw3.h>
#include <math.h>


// Reads a square image in 8-bit/color PPM format from the given file. Note: No checks on valid format are done.
void read_image(char *image_name, int pixels, double *red, double *green, double *blue)
{
    FILE *image;

    if((image = fopen(image_name, "r")))
    {
        int width, height, maxvalue;
        fscanf(image, "P6 %d %d %d ", &width, &height, &maxvalue);

        for(int row = 0; row < pixels; row++)
        for(int col = 0; col < pixels; col++)
        {
            unsigned char rgb[3];
            fread(rgb, 3, sizeof(char), image);

            red[row*pixels + col] = rgb[0];
            green[row*pixels + col] = rgb[1];
            blue[row*pixels + col] = rgb[2];
        }
        fclose(image);
    }
    else
    {
        printf("Failed to read image. '%s' could not be opened.\n", image_name);
        exit(1);
    }
}


// Writes a square image in 8-bit/color PPM format.
void write_image(char *image_name, int pixels, double *red, double *green, double *blue)
{
    FILE *image;

    if((image = fopen(image_name, "w")))
    {
        fprintf(image, "P6\n%d %d\n%d\n", pixels, pixels, 255);

        for(int row = 0; row < pixels; row++)
        for(int col = 0; col < pixels; col++)
        {
            unsigned char rgb[3];
            rgb[0] = red[row*pixels + col];
            rgb[1] = green[row*pixels + col];
            rgb[2] = blue[row*pixels + col];

            fwrite(rgb, 3, sizeof(char), image);
        }
        fclose(image);
    }
    else
    {
        printf("Failed to write image. '%s' could not be opened.\n", image_name);
        exit(1);
    }
}


int main(int argc, char **argv)
{
    const int pixels = 512, hsml = 10; // length of image in number of pixels, smoothing length h in number of pixels
    double *red, *green, *blue, *color;

    // allocate some storage for the image, and then read it
    red   = malloc(pixels * pixels * sizeof(double));
    green = malloc(pixels * pixels * sizeof(double));
    blue  = malloc(pixels * pixels * sizeof(double));

    read_image("original-image.ppm", pixels, red, green, blue);

    double total_red = 0, total_green = 0, total_blue = 0;

    for(int row = 0; row < pixels; row++)
    for(int col = 0; col < pixels; col++)
    {
        total_red += red[row*pixels + col];
        total_green += green[row*pixels + col];
        total_blue += blue[row*pixels + col];
    }
    printf("original amount of\n - red: %.25f\n - green: %.25f\n - blue: %.25f\n", total_red, total_green, total_blue);


    // Now we set up our desired smoothing kernel. We'll use complex numbers for it even though it is real.
    // first, some memory allocation
    fftw_complex *kernel_xspace = malloc(pixels * pixels * sizeof(fftw_complex));
    fftw_complex *kernel_kspace = malloc(pixels * pixels * sizeof(fftw_complex));

    // now set the values of the kernel
    double kernel_sum = 0;
    for(int row = 0; row < pixels; row++)
    for(int col = 0; col < pixels; col++)
    {
        kernel_xspace[row*pixels + col][0] = 0;  // real part
        kernel_xspace[row*pixels + col][1] = 0;  // imaginary part

        // do something sensible here to set the real part of the kernel
        double dist_from_center = sqrt(row * row + col * col);

        if (0 < dist_from_center && dist_from_center < hsml/2) {
            kernel_xspace[row*pixels + col][0] = 1 - 6 * (dist_from_center * dist_from_center)/(hsml * hsml) + 6 * (dist_from_center * dist_from_center * dist_from_center)/(hsml * hsml * hsml);
        }
        else if (hsml/2 < dist_from_center && dist_from_center < hsml) {
            kernel_xspace[row*pixels + col][0] = 2 * (1 - dist_from_center/hsml) * (1 - dist_from_center/hsml) * (1 - dist_from_center/hsml);
        }
        else {
            kernel_xspace[row*pixels + col][0] = 0;
        }
        kernel_sum += kernel_xspace[row*pixels + col][0];
    }

    for(int row = 0; row < pixels; row++)
    for(int col = 0; col < pixels; col++)
    {
        kernel_xspace[row*pixels + col][0] /= kernel_sum;
    }


    // Let's calculate the Fourier transform of the kernel
    // the FFTW3 library used here requires first a 'plan' creation, followed by the actual execution

    fftw_plan plan_kernel = fftw_plan_dft_2d (pixels, pixels, kernel_xspace, kernel_kspace, FFTW_FORWARD, FFTW_ESTIMATE);

    // now do the actual transform
    fftw_execute(plan_kernel);


    // further space allocations for image transforms
    fftw_complex *color_xspace = malloc(pixels * pixels * sizeof(fftw_complex));
    fftw_complex *color_kspace = malloc(pixels * pixels * sizeof(fftw_complex));

    // create corresponding FFT plans
    fftw_plan plan_forward =  fftw_plan_dft_2d (pixels, pixels, color_xspace, color_kspace, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_backward = fftw_plan_dft_2d (pixels, pixels, color_kspace, color_xspace, FFTW_BACKWARD, FFTW_ESTIMATE);


    // we now convolve each color channel with the kernel using FFTs
    for(int colindex = 0; colindex < 3; colindex++)
    {
        if(colindex == 0)
        color = red;
        else if(colindex == 1)
        color = green;
        else
        color = blue;

        // copy input color into complex array
        for(int row = 0; row < pixels; row++)
        for(int col = 0; col < pixels; col++)
        {
            color_xspace[row*pixels + col][0] = color[row*pixels + col];    // real part
            color_xspace[row*pixels + col][1] = 0;                      // imaginary part
        }

        // forward transform
        fftw_execute(plan_forward);

        // multiply with kernel in Fourier space
        for(int row = 0; row < pixels; row++)
        for(int col = 0; col < pixels; col++)
        {
            double stor1 = kernel_kspace[row*pixels + col][0] * color_kspace[row*pixels + col][0] - kernel_kspace[row*pixels + col][1] * color_kspace[row*pixels + col][1];


            color_kspace[row*pixels + col][1] = kernel_kspace[row*pixels + col][0] * color_kspace[row*pixels + col][1] + kernel_kspace[row*pixels + col][1] * color_kspace[row*pixels + col][0];
            color_kspace[row*pixels + col][0] = stor1;
        }

        // backward transform
        fftw_execute(plan_backward);

        // copy real value of complex result back into color array
        for(int row = 0; row < pixels; row++)
        for(int col = 0; col < pixels; col++)
        color[row*pixels + col] = color_xspace[row*pixels + col][0] / (pixels * pixels);
    }

    write_image("smoothed-image.ppm", pixels, red, green, blue);


    total_red = 0, total_green = 0, total_blue = 0;

    for(int row = 0; row < pixels; row++)
    for(int col = 0; col < pixels; col++)
    {
        total_red += red[row*pixels + col];
        total_green += green[row*pixels + col];
        total_blue += blue[row*pixels + col];
    }
    printf("smoothed amount of\n - red: %.25f\n - green: %.25f\n - blue: %.25f\n", total_red, total_green, total_blue);
    exit(0);
}
