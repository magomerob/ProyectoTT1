// $Source$
/**
 * @file global.cpp
 * @author Marcos Gómez Robres
 * @brief Esta es la implementación de la clase global.
 * @date 2025-04-23
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

 
 #include "../include/global.h"
 

Matrix eopdata;

void eop19620101(int c){
    eopdata = zeros(13,c);

    FILE *fid = fopen("eop19620101.txt", "r");
    if(fid == NULL){
        printf("Fail opening eop19620101.txt\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 1; i<=c; i++){
        fscanf(fid,"%i %d %d %i %f %f %f %f %f %f %f %f %i",&eopdata(1,i),&eopdata(2,i),&eopdata(3,i),&eopdata(4,i),&eopdata(5,i),&eopdata(6,i),&eopdata(7,i),&eopdata(8,i),&eopdata(9,i),&eopdata(10,i),&eopdata(11,i),&eopdata(12,i),&eopdata(13,i));
    }
    fclose(fid);
}