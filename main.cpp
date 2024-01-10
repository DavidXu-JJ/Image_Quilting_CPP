// Copyright 2024 Jingwei Xu. All Rights Reserved.
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <unistd.h>
#include <filesystem>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "third_party/stb_image.h"
#include "third_party/stb_image_write.h"

int get_stb_index(const int & h, const int & w, const int & width, const int & channels){
    return (width * h + w) * channels;
}

void find_patch_vertical(const unsigned char *ref, const int & ref_h, const int & ref_w,
                           const unsigned char *img, const int & height, const int & width,
                           const int & patch_size, const int & overlap_size,  const int & channels,
                           int & pick_h, int & pick_w, const float & tol, const int & max_iter){
    float min_error = 1e10;
    int min_start_h = 0;
    int min_start_w = 0;
    int offset = patch_size - overlap_size;

    bool if_finish = false;
    int iter_count = 0;
    while(!if_finish){
        int randomH = 1.0 * rand() / RAND_MAX * (height - patch_size);
        int randomW = 1.0 * rand() / RAND_MAX * (width - patch_size);
        min_start_h = randomH;
        min_start_w = randomW;
        float error = 0;
        for(int k = 0; k < patch_size; k++){
            for(int l = 0; l < overlap_size; l++){
                for(int m = 0; m < channels; m++){
                    error += pow(1.0 / 255.0 * (ref[get_stb_index(k, offset + l, ref_w, channels) + m] - img[get_stb_index(randomH + k, randomW + l, width, channels) + m]), 2);
                }
            }
        }
        error /= patch_size * overlap_size;
        if(error < min_error){
            min_error = error;
            if(error < tol){
                if_finish = true;
            }
        }
        iter_count += 1;
        if(iter_count > max_iter){
            if_finish = true;
        }
    }

    // Brute Force
    /*
    for(int i = 0; i < height - patch_size; i++){
        bool if_finish = false;
        for(int j = 0; j < width - patch_size; j++){
            float error = 0;
            for(int k = 0; k < patch_size; k++){
                for(int l = 0; l < overlap_size; l++){
                    for(int m = 0; m < channels; m++){
                        error += pow(1.0 / 255.0 * (ref[get_stb_index(k, offset + l, ref_w, channels) + m] - img[get_stb_index(i + k, j + l, width, channels) + m]), 2);
                    }
                }
            }
            error /= patch_size * overlap_size;
            if(error < min_error){
                min_error = error;
                min_start_h = i;
                min_start_w = j;
                if(error < tol){
                    if_finish = true;
                    break;
                }
            }
        }
        if(if_finish){
            break;
        }
    }
    */
    pick_h = min_start_h;
    pick_w = min_start_w;
}

std::vector<int> vertical_get_min_path(const unsigned char *ref, const int & ref_h, const int & ref_w,
                                         const unsigned char *patch_img, const int & height, const int & width,
                                         const int & patch_size, const int & overlap_size, const int & channels)
{
    int offset = patch_size - overlap_size;
    std::vector<int> min_path(patch_size);
    std::vector<std::vector<float>> dp(patch_size, std::vector<float>(overlap_size));
    for(int i = 0; i < patch_size; i++){
        for(int j = 0; j < overlap_size; j++){
            dp[i][j] = 1e10;
        }
    }
    auto get_dp = [&](const int & i, const int & j)-> float {
        if(i < 0 || i >= patch_size || j < 0 || j >= overlap_size){
            return 1e10;
        }
        return dp[i][j];
    };
    std::vector<std::vector<float>> error_vector(patch_size, std::vector<float>(overlap_size));
    for(int i = 0; i < patch_size; i++){
        for(int j = 0; j < overlap_size; j++) {
            float error = 0;
            for(int l = 0; l < channels; l++){
                error += pow(1.0 / 255.0 * (ref[get_stb_index(i, offset + j, ref_w, channels) + l] - patch_img[get_stb_index(i, j, width, channels) + l] ) , 2);
            }
            error_vector[i][j] = error;
        }
    }
    for(int i = 0; i < patch_size; i++){
        for(int j = 0; j < overlap_size; j++){
            if(i == 0){
                for(int l = 0; l < channels; l++) {
                    dp[i][j] = error_vector[i][j];
                }
            }
            else{
                for(int k = j-1; k <= j+1; k++){
                    dp[i][j] = std::min(dp[i][j], get_dp(i - 1, k) + error_vector[i][j]);
                }
            }
        }
    }
    float min_error = 1e10;
    int min_index = 0;
    for(int i = 0; i < overlap_size; i++){
        if(dp[patch_size - 1][i] < min_error){
            min_error = dp[patch_size - 1][i];
            min_index = i;
        }
    }
    min_path[patch_size - 1] = min_index;
    for(int i = patch_size - 2; i >= 0; i--){
        for(int j = min_path[i+1] - 1; j <= min_path[i+1] + 1; j++){
            if( (get_dp(i, j) + error_vector[i+1][min_path[i+1]] - get_dp(i+1, min_path[i+1])) < 1e-2 ){
                min_path[i] = j;
                break;
            }
        }
    }
    return min_path;
}

void fill_block_vertical(unsigned char * out_start_img, const int & out_h, const int & out_w,
                           const unsigned char * patch_img, const int & height, const int & width,
                           const int & patch_size, const int & overlap_size, const int & channels, const std::vector<int> & min_path,
                           const int & offset_h, const int & offset_w)
{
    for(int i = 0; i < patch_size; ++i){
        if (i + offset_h >= out_h) {
            break;
        }
        for(int j = min_path[i]; j < patch_size; ++j){
            if (j + offset_w >= out_w) {
                break;
            }
            for(int k = 0; k < channels; ++k) {
                out_start_img[get_stb_index(i, j, out_w, channels) + k] = patch_img[get_stb_index(i, j, width, channels) + k];
            }
        }
    }
}

void find_patch_horizontal(const unsigned char *ref, const int & ref_h, const int & ref_w,
                             const unsigned char *img, const int & height, const int & width,
                             const int & patch_size, const int & overlap_size,  const int & channels,
                             int & pick_h, int & pick_w, const float & tol, const int & max_iter)
{
    float min_error = 1e10;
    int min_start_h = 0;
    int min_start_w = 0;
    int offset = patch_size - overlap_size;

    bool if_finish = false;
    int iter_count = 0;
    while(!if_finish){
        int randomH = 1.0 * rand() / RAND_MAX * (height - patch_size);
        int randomW = 1.0 * rand() / RAND_MAX * (width - patch_size);
        min_start_h = randomH;
        min_start_w = randomW;
        float error = 0;
        for(int k = 0; k < patch_size; k++){
            for(int l = 0; l < overlap_size; l++){
                for(int m = 0; m < channels; m++){
                    error += pow(1.0 / 255.0 * (ref[get_stb_index(offset + l, k, ref_w, channels) + m] - img[get_stb_index(randomH + l, randomW + k, width, channels) + m]), 2);
                }
            }
        }
        error /= patch_size * overlap_size;
        if(error < min_error){
            min_error = error;
            if(error < tol){
                if_finish = true;
            }
        }
        iter_count += 1;
        if(iter_count > max_iter){
            if_finish = true;
        }
    }

    // Brute Force
    /*
    for(int i = 0; i < height - patch_size; i++){
        bool if_finish = false;
        for(int j = 0; j < width - patch_size; j++){
            float error = 0;
            for(int k = 0; k < patch_size; k++){
                for(int l = 0; l < overlap_size; l++){
                    for(int m = 0; m < channels; m++){
                        error += pow(1.0 / 255.0 * (ref[get_stb_index(offset + l, k, ref_w, channels) + m] - img[get_stb_index(i + l, j + k, width, channels) + m]), 2);
                    }
                }
            }
            if(error < min_error){
                min_error = error;
                min_start_h = i;
                min_start_w = j;
                if(error < tol){
                    if_finish = true;
                    break;
                }
            }
        }
        if(if_finish){
            break;
        }
    }
    */
    pick_h = min_start_h;
    pick_w = min_start_w;

}

std::vector<int> horizonal_get_min_path(const unsigned char *ref, const int & ref_h, const int & ref_w,
                                          const unsigned char *patch_img, const int & height, const int & width,
                                          const int & patch_size, const int & overlap_size, const int & channels)
{
    int offset = patch_size - overlap_size;
    std::vector<int> min_path(patch_size);
    std::vector<std::vector<float>> dp(overlap_size, std::vector<float>(patch_size));
    for(int i = 0; i < overlap_size; i++){
        for(int j = 0; j < patch_size; j++){
            dp[i][j] = 1e10;
        }
    }
    auto get_dp = [&](const int & i, const int & j)-> float {
        if(i < 0 || i >= overlap_size || j < 0 || j >= patch_size){
            return 1e10;
        }
        return dp[i][j];
    };
    std::vector<std::vector<float>> error_vector(overlap_size, std::vector<float>(patch_size));
    for(int i = 0; i < overlap_size; i++){
        for(int j = 0; j < patch_size; j++) {
            float error = 0;
            for(int l = 0; l < channels; l++){
                error += pow(1.0 / 255.0 * (ref[get_stb_index(offset + i, j, ref_w, channels) + l] - patch_img[get_stb_index(i, j, width, channels) + l]), 2);
            }
            error_vector[i][j] = error;
        }
    }
    for(int i = 0; i < patch_size; i++) {
        for (int j = 0; j < overlap_size; j++) {
            if (i == 0) {
                for (int l = 0; l < channels; l++) {
                    dp[j][i] = error_vector[j][i];
                }
            } else {
                for (int k = j - 1; k <= j + 1; k++) {
                    dp[j][i] = std::min(dp[j][i], get_dp(k, i-1) + error_vector[j][i]);
                }
            }
        }
    }

    float min_error = 1e10;
    int min_index = 0;
    for(int i = 0; i < overlap_size; i++){
        if(dp[i][patch_size - 1] < min_error){
            min_error = dp[i][patch_size - 1];
            min_index = i;
        }
    }
    min_path[patch_size - 1] = min_index;
    for(int i = patch_size - 2; i >= 0; i--){
        for(int j = min_path[i+1] - 1; j <= min_path[i+1] + 1; j++){
            if((get_dp(j, i) + error_vector[min_path[i+1]][i+1] - get_dp( min_path[i+1], i+1)) < 1e-2){
                min_path[i] = j;
                break;
            }
        }
    }
    return min_path;
}

void fill_block_horizontal(unsigned char * out_start_img, const int & out_h, const int & out_w,
                           const unsigned char * patch_img, const int & height, const int & width,
                           const int & patch_size, const int & overlap_size, const int & channels, const std::vector<int> & min_path,
                           const int & offset_h, const int & offset_w)
{
    for(int i = 0; i < patch_size; ++i){
        if (i + offset_w >= out_w) {
            break;
        }
        for(int j = min_path[i]; j < patch_size; ++j){
            if (j + offset_h >= out_h) {
                break;
            }
            for(int k = 0; k < channels; ++k) {
                out_start_img[get_stb_index(j, i, out_w, channels) + k] = patch_img[get_stb_index(j, i, width, channels) + k];
            }
        }
    }
}

void find_patch_both(const unsigned char *ref, const int & ref_h, const int & ref_w,
                     const unsigned char *img, const int & height, const int & width,
                     const int & patch_size, const int & overlap_size,  const int & channels,
                     int & pick_h, int & pick_w, const float & tol, const int & max_iter){
    float min_error = 1e10;
    int min_start_h = 0;
    int min_start_w = 0;
    int offset = patch_size - overlap_size;

    bool if_finish = false;
    int iter_count = 0;
    while(!if_finish){
        int randomH = 1.0 * rand() / RAND_MAX * (height - patch_size);
        int randomW = 1.0 * rand() / RAND_MAX * (width - patch_size);
        min_start_h = randomH;
        min_start_w = randomW;
        float error = 0;
        for(int k = 0; k < patch_size; k++){
            for(int l = 0; l < overlap_size; l++){
                for(int m = 0; m < channels; m++){
                    error += pow(1.0 / 255.0 * (ref[get_stb_index(k + offset, l + offset, ref_w, channels) + m] - img[get_stb_index(min_start_h + k, min_start_w + l, width, channels) + m]), 2);
                }
            }
        }
        for(int k = overlap_size; k < patch_size; k++){
            for(int l = 0; l < overlap_size; l++){
                for(int m = 0; m < channels; m++){
                    error += pow(1.0 / 255.0 * (ref[get_stb_index(l + offset, k + offset, ref_w, channels) + m] - img[get_stb_index(min_start_h + l, min_start_w + k, width, channels) + m]), 2);
                }
            }
        }
        error /= ((patch_size + offset) * overlap_size);
        if(error < min_error){
            min_error = error;
            if(error < tol){
                if_finish = true;
            }
        }
        iter_count += 1;
        if(iter_count > max_iter){
            if_finish = true;
        }
    }

    // Brute Force
    /*
    for(int i = 0; i < height - patch_size; i++){
        bool if_finish = false;
        for(int j = 0; j < width - patch_size; j++){
            float error = 0;
            for(int k = 0; k < patch_size; k++){
                for(int l = 0; l < overlap_size; l++){
                    for(int m = 0; m < channels; m++){
                        error += pow(1.0 / 255.0 * (ref[get_stb_index(k + offset, l + offset, ref_w, channels) + m] - img[get_stb_index(i + k, j + l, width, channels) + m]), 2);
                    }
                }
            }
            for(int k = overlap_size; k < patch_size; k++){
                for(int l = 0; l < overlap_size; l++){
                    for(int m = 0; m < channels; m++){
                        error += pow(1.0 / 255.0 * (ref[get_stb_index(l + offset, k + offset, ref_w, channels) + m] - img[get_stb_index(i + l, j + k, width, channels) + m]), 2);
                    }
                }
            }
            if(error < min_error){
                min_error = error;
                min_start_h = i;
                min_start_w = j;
                if(error < tol){
                    if_finish = true;
                    break;
                }
            }
        }
        if(if_finish){
            break;
        }
    }
     */
    pick_h = min_start_h;
    pick_w = min_start_w;
}

std::pair< std::vector<int>, std::vector<int> > both_get_min_path(
        const unsigned char *ref, const int & ref_h, const int & ref_w,
        const unsigned char *patch_img, const int & height, const int & width,
        const int & patch_size, const int & overlap_size, const int & channels)
{
    int offset = patch_size - overlap_size;
    std::vector<int> min_path_horizontal(patch_size);
    std::vector<std::vector<float>> dp(overlap_size, std::vector<float>(patch_size));
    for(int i = 0; i < overlap_size; i++){
        for(int j = 0; j < patch_size; j++){
            dp[i][j] = 1e10;
        }
    }
    std::function<float(const int & , const int &)> get_dp;
    get_dp = [&](const int & i, const int & j)-> float {
        if(i < 0 || i >= overlap_size || j < 0 || j >= patch_size){
            return 1e10;
        }
        return dp[i][j];
    };
    std::vector<std::vector<float>> error_vector(overlap_size, std::vector<float>(patch_size));
    for(int i = 0; i < overlap_size; i++){
        for(int j = 0; j < patch_size; j++) {
            float error = 0;
            for(int l = 0; l < channels; l++){
                error += pow(1.0 / 255.0 * (ref[get_stb_index(offset + i, offset + j, ref_w, channels) + l] - patch_img[get_stb_index(i, j, width, channels) + l]), 2);
            }
            error_vector[i][j] = error;
        }
    }
    for(int i = 0; i < patch_size; i++) {
        for (int j = 0; j < overlap_size; j++) {
            if (i == 0) {
                for (int l = 0; l < channels; l++) {
                    dp[j][i] = error_vector[j][i];
                }
            } else {
                for (int k = j - 1; k <= j + 1; k++) {
                    dp[j][i] = std::min(dp[j][i], get_dp(k, i-1) + error_vector[j][i]);
                }
            }
        }
    }

    float min_error = 1e10;
    int min_index = 0;
    for(int i = 0; i < overlap_size; i++){
        if(dp[i][patch_size - 1] < min_error){
            min_error = dp[i][patch_size - 1];
            min_index = i;
        }
    }
    min_path_horizontal[patch_size - 1] = min_index;
    for(int i = patch_size - 2; i >= 0; i--){
        for(int j = min_path_horizontal[i+1] - 1; j <= min_path_horizontal[i+1] + 1; j++){
            if((get_dp(j, i) + error_vector[min_path_horizontal[i+1]][i+1] - get_dp( min_path_horizontal[i+1], i+1)) < 1e-2){
                min_path_horizontal[i] = j;
                break;
            }
        }
    }

    std::vector<int> min_path_vertical(patch_size);
    dp.clear();
    dp.assign(patch_size, std::vector<float>(overlap_size));
    for(int i = 0; i < patch_size; i++){
        for(int j = 0; j < overlap_size; j++){
            dp[i][j] = 1e10;
        }
    }
    get_dp = [&](const int & i, const int & j)-> float {
        if(i < 0 || i >= patch_size || j < 0 || j >= overlap_size){
            return 1e10;
        }
        return dp[i][j];
    };
    error_vector.clear();
    error_vector.assign(patch_size, std::vector<float>(overlap_size));
    for(int i = 0; i < patch_size; i++){
        for(int j = 0; j < overlap_size; j++) {
            float error = 0;
            for(int l = 0; l < channels; l++){
                error += pow(1.0 / 255.0 * (ref[get_stb_index(offset + i, offset + j, ref_w, channels) + l] - patch_img[get_stb_index(i, j, width, channels) + l] ) , 2);
            }
            error_vector[i][j] = error;
        }
    }
    for(int i = 0; i < patch_size; i++){
        for(int j = 0; j < overlap_size; j++){
            if(i == 0){
                for(int l = 0; l < channels; l++) {
                    dp[i][j] = error_vector[i][j];
                }
            }
            else{
                for(int k = j-1; k <= j+1; k++){
                    dp[i][j] = std::min(dp[i][j], get_dp(i - 1, k) + error_vector[i][j]);
                }
            }
        }
    }
    min_error = 1e10;
    min_index = 0;
    for(int i = 0; i < overlap_size; i++){
        if(dp[patch_size - 1][i] < min_error){
            min_error = dp[patch_size - 1][i];
            min_index = i;
        }
    }
    min_path_vertical[patch_size - 1] = min_index;
    for(int i = patch_size - 2; i >= 0; i--){
        for(int j = min_path_vertical[i+1] - 1; j <= min_path_vertical[i+1] + 1; j++){
            if( (get_dp(i, j) + error_vector[i+1][min_path_vertical[i+1]] - get_dp(i+1, min_path_vertical[i+1])) < 1e-2 ){
                min_path_vertical[i] = j;
                break;
            }
        }
    }
    return std::make_pair(min_path_horizontal, min_path_vertical);
}

void fill_block_both(unsigned char * out_start_img, const int & out_h, const int & out_w,
                     const unsigned char * patch_img, const int & height, const int & width,
                     const int & patch_size, const int & overlap_size, const int & channels,
                     const std::vector<int> & min_path_vertical, const std::vector<int> & min_path_horizontal,
                     const int & offset_h, const int & offset_w)
{
    for(int i = 0; i < patch_size; ++i){
        if (i + offset_w >= out_w) {
            break;
        }
        for(int j = min_path_horizontal[i]; j < patch_size; ++j){
            if (i < min_path_vertical[j]) {
                continue;
            }
            if (j + offset_h >= out_h) {
                break;
            }
            for(int k = 0; k < channels; ++k) {
                out_start_img[get_stb_index(j, i, out_w, channels) + k] = patch_img[get_stb_index(j, i, width, channels) + k];
            }
        }
    }
}



void image_quilting(
        const std::string & input_path, const std::string & output_path,
        const int & scale, const int & patch_size, const float & overlap, const float & tol, const int & n_iter
){
    // Read image
    int width, height, channels;
    unsigned char *img = stbi_load(input_path.c_str(), &width, &height, &channels, 0);
    if (img == NULL) {
        std::cout << "Error in loading the image" << std::endl;
        exit(1);
    }
    std::cout << "Image loaded successfully" << std::endl;

    int out_h = scale * height;
    int out_w = scale * width;
    int overlap_size = patch_size * overlap;

    int nH = int(ceil((out_h - patch_size) * 1.0 / (patch_size - overlap_size)));
    int nW = int(ceil((out_w - patch_size) * 1.0 / (patch_size - overlap_size)));

    out_h = patch_size + nH * (patch_size - overlap_size);
    out_w = patch_size + nW * (patch_size - overlap_size);

    unsigned char *img_out = new unsigned char[out_h * out_w * channels];
    for(int i = 0; i < out_h; ++i){
        for(int j = 0; j < out_w; ++j){
            for(int k = 0; k < channels; ++k){
                img_out[get_stb_index(i, j, out_w, channels) + k] = 0;
            }
        }
    }
    std::srand(std::time(nullptr));

    int randomH = 1.0 * rand() / RAND_MAX * (height - patch_size);
    int randomW = 1.0 * rand() / RAND_MAX * (width - patch_size);

    for (int i = 0; i < patch_size; i++) {
        for (int j = 0; j < patch_size; j++) {
            for (int k = 0; k < channels; k++) {
                img_out[get_stb_index(i, j, out_w, channels) + k] = img[get_stb_index(randomH + i, randomW + j, width, channels) + k];
            }
        }
    }

    int step_size = patch_size - overlap_size;

    for(int start_w = 0; start_w < out_w; start_w += step_size){
        int pick_h = 0, pick_w = 0;
        int current_w = std::min(start_w, out_w - patch_size - 1);

        unsigned char * ref_img = img_out + get_stb_index(0, current_w, out_w, channels);

        find_patch_vertical(ref_img, out_h, out_w,
                              img, height, width,
                              patch_size, overlap_size, channels,
                              pick_h, pick_w, tol, n_iter);


        unsigned char * patch_img = img + get_stb_index(pick_h, pick_w, width, channels);

        std::vector<int> min_path = vertical_get_min_path(ref_img, out_h, out_w,
                                                            patch_img, height, width,
                                                            patch_size, overlap_size, channels);

        unsigned char * out_start_img = img_out + get_stb_index(0, current_w + step_size + 1, out_w, channels);

        fill_block_vertical(out_start_img, out_h, out_w,
                            patch_img, height, width,
                            patch_size, overlap_size, channels, min_path, 0, current_w + step_size + 1);

    }

    for(int start_h = 0; start_h < out_h; start_h += step_size){
        int pick_h = 0, pick_w = 0;
        int current_h = std::min(start_h, out_h - patch_size - 1);

        unsigned char * ref_img = img_out + get_stb_index(current_h, 0, out_w, channels);

        find_patch_horizontal(ref_img, out_h, out_w,
                            img, height, width,
                            patch_size, overlap_size, channels,
                            pick_h, pick_w, tol, n_iter);

        unsigned char * patch_img = img + get_stb_index(pick_h, pick_w, width, channels);

        std::vector<int> min_path = horizonal_get_min_path(ref_img, out_h, out_w,
                                                          patch_img, height, width,
                                                          patch_size, overlap_size, channels);


        unsigned char * out_start_img = img_out + get_stb_index(current_h + step_size + 1, 0 , out_w, channels);

        fill_block_horizontal(out_start_img, out_h, out_w,
                              patch_img, height, width,
                              patch_size, overlap_size, channels, min_path, current_h + step_size + 1, 0);

    }

    for(int start_h = 0; start_h < out_h; start_h += step_size){
        for(int start_w = 0; start_w < out_w; start_w += step_size){
            int pick_h = 0, pick_w = 0;
            int current_h = std::min(start_h, out_h - patch_size - 1 - step_size);
            int current_w = std::min(start_w, out_w - patch_size - 1 - step_size);

            unsigned char * ref_img = img_out + get_stb_index(current_h, current_w, out_w, channels);

            find_patch_both(ref_img, out_h, out_w,
                            img, height, width,
                            patch_size, overlap_size, channels,
                            pick_h, pick_w, tol, n_iter);

            unsigned char * patch_img = img + get_stb_index(pick_h, pick_w, width, channels);

            auto both_path = both_get_min_path(ref_img, out_h, out_w,
                                               patch_img, height, width,
                                               patch_size, overlap_size, channels);

            unsigned char * out_start_img = img_out + get_stb_index(current_h + step_size + 1, current_w + step_size + 1, out_w, channels);

            fill_block_both(out_start_img, out_h, out_w,
                            patch_img, height, width,
                            patch_size, overlap_size, channels, both_path.second, both_path.first,
                            current_h + step_size + 1, current_w + step_size + 1);

        }
    }

    // Write image
    stbi_write_png(output_path.c_str(), out_w, out_h, channels, img_out, out_w * channels);
    std::cout << "Image written successfully" << std::endl;
    stbi_image_free(img);
    delete [] img_out;
}


int main(int argc, char** argv) {
    std::string input_path = "input/example.png";
    std::string output_path = "output/example.png";
    std::filesystem::path dir = std::filesystem::current_path();
    while (dir != std::filesystem::current_path().root_path() && !exists(dir / input_path))
    {
        dir = dir.parent_path();
    }
    if (dir == std::filesystem::current_path().root_path()){
        std::cout << "Something wrong happens. The image file is not found.\n";
        return 1;
    }
    input_path = dir / input_path;
    output_path = dir / output_path;

    int block = 50;
    float overlap = 1.0/6;
    int scale = 2;
    float tol = 0.01;

    int opt;
    while ((opt = getopt(argc, argv, "i:o:")) != -1) {
        switch (opt) {
            case 'i':
                input_path = optarg;
                break;
            case 'o':
                output_path = optarg;
                break;
            default:
                std::cerr << "Usage: " << argv[0] << " [-i input_path] [-o output_path]\n";
                return 1;
        }
    }

    image_quilting(input_path, output_path, scale, block, overlap, tol, 256);

    std::cout << input_path << std::endl;

    return 0;
}
