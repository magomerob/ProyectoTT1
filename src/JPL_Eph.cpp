// $Source$
/**
 * @file JPL_Eph.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de JPL_Eph.
 * @version 0.1
 * @date 2025-05-07
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

 #include "../include/JPL_Eph.h"
Matrix rango(int i, int f, int s){
    int len = (f-i)/s;
    Matrix m(len+1);
    for (size_t j = 0; j <= len; j++)
    {
        m(j+1) = (i+(s*j*1.0));
    }
    return m;
}

tuple<Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix> JPL_Eph(double Mjd_TDB)
{

    Matrix temp, PCtemp, Cy_Earth, Cz_Earth, Cx_Earth, Cx, Cy, Cz, Cy_Moon, Cz_Moon, Cx_Moon, Nutations, Cx_Sun, Cz_Sun, Cy_Sun,
    Cx_Mercury, Cy_Mercury, Cz_Mercury, Cx_Venus, Cy_Venus, Cz_Venus, Cx_Mars, Cy_Mars, Cz_Mars, Cx_Jupiter, Cy_Jupiter, Cz_Jupiter, 
    Cx_Saturn, Cy_Saturn, Cz_Saturn, Cx_Uranus, Cy_Uranus, Cz_Uranus, Cx_Neptune,Cy_Neptune,Cz_Neptune, Cx_Pluto,Cy_Pluto,Cz_Pluto, v1, v2, v3;
    double JD = Mjd_TDB + 2400000.5;
    double i;
    for(i=1.0; i < PC.fil; i++) {
        if (PC(i, 1)<=JD && JD<=PC(i, 2)){
            break;
        } 
    }
PCtemp = PC.extract_row(i);

double t1 = PCtemp(1)-2400000.5;

double dt = Mjd_TDB - t1;

temp = rango(231,270,13);

Cx_Earth = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Earth = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Earth = extract_vector(PCtemp,temp(3),temp(4)-1);
temp = temp+39;
Cx = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz = extract_vector(PCtemp,temp(3),temp(4)-1);
Cx_Earth = union_vector(Cx_Earth,Cx);
Cy_Earth = union_vector(Cy_Earth,Cy);
Cz_Earth = union_vector(Cz_Earth,Cz);
double j, Mjd0;
if (0<=dt && dt<=16){
    j=0;
    Mjd0 = t1;
}else if(16<dt && dt<=32){
    j=1;
    Mjd0 = t1+16*j;
}

v1 = extract_vector(Cx_Earth,13*j+1,13*j+13);
v2 = extract_vector(Cy_Earth,13*j+1,13*j+13);
v3 = (extract_vector(Cz_Earth,13*j+1,13*j+13));

Matrix r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, v1,v2, v3)*1e3;
//cout<< "r_Earth: " << r_Earth << endl;
temp = rango(441,480,13);
Cx_Moon = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Moon = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Moon = extract_vector(PCtemp,temp(3),temp(4)-1);
for (i=1; i<=7; i++){
    temp = temp+39;
    Cx = extract_vector(PCtemp,temp(1),temp(2)-1);
    Cy = extract_vector(PCtemp,temp(2),temp(3)-1);
    Cz = extract_vector(PCtemp,temp(3),temp(4)-1);   
    Cx_Moon = union_vector(Cx_Moon,Cx);
    Cy_Moon = union_vector(Cy_Moon,Cy);
    Cz_Moon = union_vector(Cz_Moon,Cz);    
}
if (0<=dt && dt<=4){
    j=0;
    Mjd0 = t1;
}else if(4<dt && dt<=8){
    j=1;
    Mjd0 = t1+4*j;
}else if(8<dt && dt<=12){
    j=2;
    Mjd0 = t1+4*j;
}else if(12<dt && dt<=16){
    j=3;
    Mjd0 = t1+4*j;
}else if(16<dt && dt<=20){
    j=4;
    Mjd0 = t1+4*j;
}else if(20<dt && dt<=24){
    j=5;
    Mjd0 = t1+4*j;
}else if(24<dt && dt<=28){
    j=6;
    Mjd0 = t1+4*j;
}else if(28<dt && dt<=32){
    j=7;
    Mjd0 = t1+4*j;
}
v1 = extract_vector(Cx_Moon,13*j+1,13*j+13);
v2 = extract_vector(Cy_Moon,13*j+1,13*j+13);
v3 = (extract_vector(Cz_Moon,13*j+1,13*j+13));
Matrix r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, v1,v2,v3)*1e3;
//cout<< "r_Moon: " << r_Moon << endl;
temp = rango(753,786,11);
Cx_Sun = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Sun = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Sun = extract_vector(PCtemp,temp(3),temp(4)-1);
temp = temp+33;
Cx = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz = extract_vector(PCtemp,temp(3),temp(4)-1);   
Cx_Sun = union_vector(Cx_Sun,Cx);
Cy_Sun = union_vector(Cy_Sun,Cy);
Cz_Sun = union_vector(Cz_Sun,Cz);
if (0<=dt && dt<=16){
    j=0;
    Mjd0 = t1;
}else if(16<dt && dt<=32){
    j=1;
    Mjd0 = t1+16*j;
}
v1 = extract_vector(Cx_Sun,11*j+1,11*j+11);
v2 = extract_vector(Cy_Sun,11*j+1,11*j+11);
v3 = (extract_vector(Cz_Sun,11*j+1,11*j+11));
Matrix r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16, v1,v2,v3)*1e3;
//cout<< "r_Sun: " << r_Sun << endl;
temp = rango(3,45,14);
Cx_Mercury = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Mercury = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Mercury = extract_vector(PCtemp,temp(3),temp(4)-1);
for(int i=1; i<=3; i++){
    temp = temp+42;
    Cx = extract_vector(PCtemp,temp(1),temp(2)-1);
    Cy = extract_vector(PCtemp,temp(2),temp(3)-1);
    Cz = extract_vector(PCtemp,temp(3),temp(4)-1);
    Cx_Mercury = union_vector(Cx_Mercury,Cx);
    Cy_Mercury = union_vector(Cy_Mercury,Cy);
    Cz_Mercury = union_vector(Cz_Mercury,Cz);
}
if (0<=dt && dt<=8){
    j=0;
    Mjd0 = t1;
}else if(8<dt && dt<=16){
    j=1;
    Mjd0 = t1+8*j;
}else if (16<dt && dt<=24){
    j=2;
    Mjd0 = t1+8*j;
}else if(24<dt && dt<=32){
    j=3;
    Mjd0 = t1+8*j;
}
v1 = extract_vector(Cx_Mercury,14*j+1,14*j+14);
v2 = extract_vector(Cy_Mercury,14*j+1,14*j+14);
v3 = (extract_vector(Cz_Mercury,14*j+1,14*j+14));
Matrix r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, v1,v2,v3)*1e3;
//cout<< "r_Mercury: " << r_Mercury << endl;
temp = rango(171,201,10);
Cx_Venus = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Venus = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Venus = extract_vector(PCtemp,temp(3),temp(4)-1);
temp = temp+30;
Cx = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz = extract_vector(PCtemp,temp(3),temp(4)-1);
Cx_Venus = union_vector(Cx_Venus,Cx);
Cy_Venus = union_vector(Cy_Venus,Cy);
Cz_Venus = union_vector(Cz_Venus,Cz);
if (0<=dt && dt<=16){
    j=0;
    Mjd0 = t1;
}else if(16<dt && dt<=32){
    j=1;
    Mjd0 = t1+16*j;
}
v1 = extract_vector(Cx_Venus,10*j+1,10*j+10);
v2 = extract_vector(Cy_Venus,10*j+1,10*j+10);
v3 = (extract_vector(Cz_Venus,10*j+1,10*j+10));
Matrix r_Venus = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, v1,v2,v3)*1e3;
//cout<< "r_Venus: " << r_Venus << endl;
temp = rango(309,342,11);
Cx_Mars = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Mars = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Mars = extract_vector(PCtemp,temp(3),temp(4)-1);
j=0;
Mjd0 = t1;
v1 = extract_vector(Cx_Mars,11*j+1,11*j+11);
v2 = extract_vector(Cy_Mars,11*j+1,11*j+11);
v3 = (extract_vector(Cz_Mars,11*j+1,11*j+11));
Matrix r_Mars = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, v1,v2,v3)*1e3;
//cout<< "r_Mars: " << r_Mars << endl;
temp = rango(342,366,8);
Cx_Jupiter = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Jupiter = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Jupiter = extract_vector(PCtemp,temp(3),temp(4)-1);
j=0;
Mjd0 = t1;
v1 = extract_vector(Cx_Jupiter,8*j+1,8*j+8);
v2 = extract_vector(Cy_Jupiter,8*j+1,8*j+8);
v3 = (extract_vector(Cz_Jupiter,8*j+1,8*j+8));
Matrix r_Jupiter = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32, v1,v2,v3)*1e3;
//cout<< "r_Jupiter: " << r_Jupiter << endl;
temp = rango(366,387,7);
Cx_Saturn = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Saturn = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Saturn = extract_vector(PCtemp,temp(3),temp(4)-1);
j=0;
Mjd0 = t1;
v1 = extract_vector(Cx_Saturn,7*j+1,7*j+7);
v2 = extract_vector(Cy_Saturn,7*j+1,7*j+7);
v3 = (extract_vector(Cz_Saturn,7*j+1,7*j+7));

Matrix r_Saturn = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32, v1,v2,v3)*1e3;
//cout<< "r_Saturn: " << r_Saturn << endl;    
temp = rango(387,405,6);
Cx_Uranus = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Uranus = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Uranus = extract_vector(PCtemp,temp(3),temp(4)-1);
j=0;
Mjd0 = t1;
v1 = extract_vector(Cx_Uranus,6*j+1,6*j+6);
v2 = extract_vector(Cy_Uranus,6*j+1,6*j+6);
v3 = (extract_vector(Cz_Uranus,6*j+1,6*j+6));
Matrix r_Uranus = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, v1,v2,v3)*1e3;
//cout<< "r_Uranus: " << r_Uranus << endl;
temp = rango(405,423,6);
Cx_Neptune = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Neptune = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Neptune = extract_vector(PCtemp,temp(3),temp(4)-1);
j=0;
Mjd0 = t1;
v1 = extract_vector(Cx_Neptune,6*j+1,6*j+6);
v2 = extract_vector(Cy_Neptune,6*j+1,6*j+6);
v3 = (extract_vector(Cz_Neptune,6*j+1,6*j+6));
Matrix r_Neptune = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, v1,v2,v3)*1e3;
//cout<< "r_Neptune: " << r_Neptune << endl;
temp = rango(423,441,6);
Cx_Pluto = extract_vector(PCtemp,temp(1),temp(2)-1);
Cy_Pluto = extract_vector(PCtemp,temp(2),temp(3)-1);
Cz_Pluto = extract_vector(PCtemp,temp(3),temp(4)-1);
j=0;
Mjd0 = t1;
v1 = extract_vector(Cx_Pluto,6*j+1,6*j+6);
v2 = extract_vector(Cy_Pluto,6*j+1,6*j+6);
v3 = (extract_vector(Cz_Pluto,6*j+1,6*j+6));

Matrix r_Pluto = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, v1,v2,v3)*1e3;
//cout<< "r_Pluto: " << r_Pluto << endl;
/*
temp = rango(819,839,10);
Cx_Nutations = PCtemp(temp(1):temp(2)-1);
Cy_Nutations = PCtemp(temp(2):temp(3)-1);
for i=1:3
    temp = temp+20;
    Cx = PCtemp(temp(1):temp(2)-1);
    Cy = PCtemp(temp(2):temp(3)-1);
    Cx_Nutations = [Cx_Nutations,Cx];
    Cy_Nutations = [Cy_Nutations,Cy];
end
if (0<=dt && dt<=8)
    j=0;
    Mjd0 = t1;
elseif(8<dt && dt<=16)
    j=1;
    Mjd0 = t1+8*j;
elseif (16<dt && dt<=24)
    j=2;
    Mjd0 = t1+8*j;
elseif(24<dt && dt<=32)
    j=3;
    Mjd0 = t1+8*j;
end

Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Nutations(10*j+1:10*j+10),...
                   Cy_Nutations(10*j+1:10*j+10),zeros(10,1))';

temp = (899:10:929);
Cx_Librations = PCtemp(temp(1):temp(2)-1);
Cy_Librations = PCtemp(temp(2):temp(3)-1);
Cz_Librations = PCtemp(temp(3):temp(4)-1);
for i=1:3
    temp = temp+30;
    Cx = PCtemp(temp(1):temp(2)-1);
    Cy = PCtemp(temp(2):temp(3)-1);
    Cz = PCtemp(temp(3):temp(4)-1);
    Cx_Librations = [Cx_Librations,Cx];
    Cy_Librations = [Cy_Librations,Cy];
    Cz_Librations = [Cz_Librations,Cz];    
end
if (0<=dt && dt<=8)
    j=0;
    Mjd0 = t1;
elseif(8<dt && dt<=16)
    j=1;
    Mjd0 = t1+8*j;
elseif (16<dt && dt<=24)
    j=2;
    Mjd0 = t1+8*j;
elseif(24<dt && dt<=32)
    j=3;
    Mjd0 = t1+8*j;
end
Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Librations(10*j+1:10*j+10),...
                    Cy_Librations(10*j+1:10*j+10), Cz_Librations(10*j+1:10*j+10))';
*/

double EMRAT = 81.30056907419062;
double EMRAT1 = 1/(1+EMRAT);

r_Earth = r_Earth-(r_Moon*EMRAT1);
r_Mercury = (r_Earth*-1)+r_Mercury;
r_Venus = (r_Earth*-1)+r_Venus;
r_Mars = (r_Earth*-1)+r_Mars;
r_Jupiter = (r_Earth*-1)+r_Jupiter;
r_Saturn = (r_Earth*-1)+r_Saturn;
r_Uranus = (r_Earth*-1)+r_Uranus;
r_Neptune = (r_Earth*-1)+r_Neptune;
r_Pluto = (r_Earth*-1)+r_Pluto;
r_Sun = (r_Earth*-1)+r_Sun;
return make_tuple(r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, r_Neptune,r_Pluto,r_Moon,r_Sun);
}

 