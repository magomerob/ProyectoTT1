// $Source$
/**
 * @file accel.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de la aceleración.
 * @version 0.1
 * @date 2025-04-30
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/accel.h"

Matrix accel(double x, Matrix y)
{
    auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = iers(eopdata,AuxParam.Mjd_UTC + x/86400,'l');
    auto [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
    long double Mjd_UT1 = AuxParam.Mjd_UTC + x/86400 + UT1_UTC/86400;
    long double Mjd_TT = AuxParam.Mjd_UTC + x/86400 + TT_UTC/86400;
    
    Matrix P = precMatrix(SAT_Const::MJD_J2000,Mjd_TT);
    
    Matrix N = nutMatrix(Mjd_TT);
    
    Matrix T = N * P;
    
    Matrix E = poleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    long double MJD_TDB = Mjd_TDB(Mjd_TT);
    auto [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, r_Neptune,r_Pluto,r_Moon,r_Sun] = JPL_Eph(MJD_TDB);

    // Acceleration due to harmonic gravity field
    Matrix v = transpose(extract_vector(transpose(y),1,3));
    
    Matrix a = accelHarmonic(v, E, AuxParam.n, AuxParam.m);
    //cout<<"a:"<<endl<<a<<endl;
    // Luni-solar perturbations
    Matrix vt = transpose(v);
    if (AuxParam.sun){
        a = a + transpose(AccelPointMass(vt,r_Sun,SAT_Const::GM_Sun));
    }

    if (AuxParam.moon){
        a = a + transpose(AccelPointMass(vt,r_Moon,SAT_Const::GM_Moon));
    }

    // Planetary perturbations
    
    if (AuxParam.planets){
        a = a + transpose(AccelPointMass(vt,r_Mercury,SAT_Const::GM_Mercury));
        a = a + transpose(AccelPointMass(vt,r_Venus,SAT_Const::GM_Venus));
        a = a + transpose(AccelPointMass(vt,r_Mars,SAT_Const::GM_Mars));
        a = a + transpose(AccelPointMass(vt,r_Jupiter,SAT_Const::GM_Jupiter));
        a = a + transpose(AccelPointMass(vt,r_Saturn,SAT_Const::GM_Saturn));
        a = a + transpose(AccelPointMass(vt,r_Uranus,SAT_Const::GM_Uranus));
        a = a + transpose(AccelPointMass(vt,r_Neptune,SAT_Const::GM_Neptune));
        a = a + transpose(AccelPointMass(vt,r_Pluto,SAT_Const::GM_Pluto));
    }

    Matrix v2 = extract_vector(transpose(y),4,6);

    return transpose(union_vector(v2,transpose(a)));
    //dY = [Y(4:6);a];

}