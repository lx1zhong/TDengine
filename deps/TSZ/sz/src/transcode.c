/**
 * @file transcode.c
 * @author Yu Zhong (lx1zhong@qq.com)
 * @brief 
 * @date 2022-08-29
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "transcode.h"
#include "fse.h"
#include "bitstream.h"
#include <stdlib.h>


// int comp(const void *a, const void *b) {
//     return *(int *)a - *(int *)b;
// }

/**
 * @brief transform type array to FseCode & tranCodeBits
 * [type] ---minus md---> [factor] ---transcode---> [tp_code]+[bitstream of diff]
 * 
 */
void encode_with_fse(int *type, size_t dataSeriesLength, unsigned int intervals, 
                    unsigned char **FseCode, size_t *FseCode_size, 
                    unsigned char **transCodeBits, size_t *transCodeBits_size) {
    // int max_type=0, min_type=0;
    // int diff;

    // int *type_ = (int*)malloc(dataSeriesLength*sizeof(int));
    // memcpy(type_, type, dataSeriesLength*sizeof(int));

    // qsort(type_, dataSeriesLength, sizeof(int), comp);
    // int first=0;
    // while(type_[first]==0) first++;
    // printf("【type】max=%d, mid=%d, intervals/2=%d\n", type_[dataSeriesLength-1], type_[(dataSeriesLength-first)/2 + first], intervals/2);

    // transcoding results of type array
    uint8_t *tp_code = (uint8_t *)malloc(dataSeriesLength);
    // uint8_t tp_code[dataSeriesLength];           // 数值太大时，初始化数组会失败(why?)
    int factor;
    int nbits;

    BIT_CStream_t transCodeStream;
    int dstCapacity = dataSeriesLength * sizeof(int);
    void *tmp = malloc(dstCapacity);
    BIT_initCStream(&transCodeStream, tmp, dstCapacity);

    // transcoding
    // printf("intervals=%u\n", intervals);
    int md = intervals/2;
    for (int i=0; i<dataSeriesLength; i++) {
        // if (type[i] > max_type)
        // 	max_type = type[i];
        // if (type[i] < min_type)
        // 	min_type = type[i];
        if (type[i] == 0) {
            // unpredictable data
            tp_code[i] = 67;
            nbits = 0;
        }
        else {
            factor = type[i] - md;
            if(factor >= 32768 || factor <= -32768) {
                printf("wrong! type=%d\n", type[i]);
                exit(1); 
            }
            tp_code[i] = (uint8_t)Int2code(factor);
            nbits = code2int[tp_code[i]][1];
            if (factor >= 0) {
                BIT_addBits(&transCodeStream, factor, nbits);
                // diff = (factor & BIT_mask[nbits]);
            }
            else {
                BIT_addBits(&transCodeStream, -factor, nbits);
                // diff = ((-factor) & BIT_mask[nbits]);
            }
            // printf(" %d: factor=%d, nbits=%d, tp_code=%d, diff=%d\n", i, factor, nbits, tp_code[i], diff);
            BIT_flushBits(&transCodeStream);
        }
        
    }
    
    // fse encoding
    (*FseCode) = (unsigned char*)malloc(dataSeriesLength);
    size_t fse_size = FSE_compress((*FseCode), dataSeriesLength, tp_code, dataSeriesLength);
    if (FSE_isError(fse_size)) {
        printf("encode:FSE_isError!\n");
        exit(1);
    }
    printf("fse_size=%lu, ", fse_size);
    (*FseCode_size) = fse_size;

    size_t const streamSize = BIT_closeCStream(&transCodeStream);
    printf("streamSize=%lu\n", streamSize);
    (*transCodeBits_size) = streamSize;
    if (streamSize == 0) {
        printf("too small!\n");   /* not enough space */
        exit(1);
    }
    (*transCodeBits) = malloc(streamSize);
    memcpy((*transCodeBits), tmp, streamSize);
    free(tmp);

    // FILE *f0 = fopen("/home/lxzhong/tmp/type_array.bin","wb");
    // fwrite(type, sizeof(int), dataSeriesLength, f0);
    // fclose(f0);

    // FILE *f1 = fopen("/home/lxzhong/tmp/tp_code.txt","w");
    // for (int j=0; j<dataSeriesLength; j++) {
    // 	fprintf(f1, "%d ", (int)tp_code[j]);
    // }
    // fclose(f1);

    // FILE *f2 = fopen("/home/lxzhong/tmp/fse.bin","wb");
    // fwrite((*FseCode), 1, fse_size, f2);
    // fclose(f2);

    // FILE *f3 = fopen("/home/lxzhong/tmp/transCodeBits.bin","wb");
    // fwrite((*transCodeBits), 1, streamSize, f3);
    // fclose(f3);

    // int* type2 = (int*)malloc(dataSeriesLength*sizeof(int));
    // decode_with_fse((*this), dataSeriesLength, type2);
    // printf("max_type=%d, min_type=%d\n", max_type, min_type);

    free(tp_code);
}

/**
 * @brief transform FseCode & tranCodeBits to type array 
 * 
 */
void decode_with_fse(int *type, size_t dataSeriesLength, unsigned int intervals, 
                    unsigned char *FseCode, size_t FseCode_size, 
                    unsigned char *transCodeBits, size_t transCodeBits_size) {

    uint8_t *tp_code = (uint8_t *)malloc(dataSeriesLength);

    size_t fse_size = FSE_decompress(tp_code, dataSeriesLength, FseCode, FseCode_size);
    if (FSE_isError(fse_size)) {
        printf("decode:FSE_isError!\n");
        exit(1);
    }
    if (fse_size != dataSeriesLength) {
        printf("fse_size(%lu) != dataSeriesLength(%lu)!\n", fse_size, dataSeriesLength);
        exit(1);
    }

    BIT_DStream_t transCodeStream;
    size_t stream_size = BIT_initDStream(&transCodeStream, transCodeBits, transCodeBits_size);
    if (stream_size == 0) {
        printf("transcode stream empty!\n");
        exit(1);
    }

    int md = intervals / 2;
    // printf("intervals = %u\n", intervals);

    int factor, nbits, base;
    size_t diff;
    for (int i=dataSeriesLength-1; i>=0; i--) {
        if (tp_code[i] == 67) {
            type[i] = 0;
            continue;
        }
        base = code2int[tp_code[i]][0];
        nbits = code2int[tp_code[i]][1];

        if (nbits >= 1)
            diff = BIT_readBitsFast(&transCodeStream, nbits);
        else
            diff =0;
        BIT_reloadDStream(&transCodeStream);
        
        if (base >= 0) {
            factor = base + diff;
        }
        else if (base < 0) {
            factor = base - diff;
        }
        type[i] = factor + md;
		// printf(" %d:factor=%d, base=%d, nbits=%d, tp_code=%d, diff=%d\n", i, factor, base, nbits, tp_code[i], diff);
    }

    
    // FILE *f0 = fopen("/home/lxzhong/tmp/type_array2.bin","wb");
    // fwrite(type, sizeof(int), dataSeriesLength, f0);
    // fclose(f0);

}

