/*
 *   FAC - Flexible Atomic Code
 *   Copyright (C) 2001-2015 Ming Feng Gu
 * 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 * 
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 * 
 *   You should have received a copy of the GNU General Public License
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "nucleus.h"
#include "cf77.h"
#include "coulomb.h"

static char *rcsid="$Id: nucleus.c,v 1.14 2005/01/17 05:39:41 mfgu Exp $";
#if __GNUC__ == 2
#define USE(var) static void * use_##var = (&use_##var, (void *) &var) 
USE (rcsid);
#endif

static NUCLEUS atom;
static char _ename[N_ELEMENTS][3] = 
{"H", "He", "Li", "Be", "B", "C", "N", "O", "F",
 "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", 
 "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", 
 "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
 "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", 
 "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
 "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", 
 "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
 "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
 "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
 "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
 "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs",
 "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og", "Ux", "Uy"};

static double _emass[N_ELEMENTS] = 
{1.008, 4.003, 6.94, 9.012, 10.81, 12.01, 14.01, 16.00, 19.00,
 20.18, 22.99, 24.31, 26.98, 28.09, 30.97, 32.06, 35.45, 39.95, 39.10, 
 40.08, 44.96, 47.88, 50.94, 52.00, 54.94, 55.85, 58.93, 58.69, 63.55,
 65.39, 69.72, 72.64, 74.92, 78.96, 79.90, 83.79, 85.47,
 87.62, 88.91, 91.22, 92.91, 95.96, 98.00, 101.1, 102.9,
 106.4, 107.9, 112.4, 114.8, 118.7, 121.8, 127.6, 126.9,
 131.3, 132.9, 137.3, 138.9, 140.1, 140.9, 144.2, 145.0,
 150.4, 152.0, 157.2, 158.9, 162.500001, 164.9, 167.3,
 168.9, 173.0, 175.0, 178.5, 180.9, 183.9, 186.2, 190.2,
 192.2, 195.1, 197.0, 200.5, 204.38, 207.2, 209.0, 
 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0,
 231.0, 238.0, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 
 252.0, 257.0, 258.0, 259.0, 262.0, 267.0, 268.0, 269.0, 270.0,
 277.0, 278., 281.0, 282.0, 285.0, 286.0, 289.0, 289.0, 293.0, 294.0, 294.0, 296.0, 296.0};

static double _arrms[N_ELEMENTS] =
  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
   3.057, 3.061, 3.122, 3.189, 3.261, 3.365, 3.427, 3.435,
   3.476, 3.544, 3.591, 3.599, 3.642, 3.706, 3.737, 3.788,
   3.775, 3.882, 3.929, 3.997, 4.074, 4.097, 4.140, 4.163,
   4.188, 4.203, 4.220, 4.242, 4.270, 4.324, 4.409, 4.424,
   4.482, 4.494, 4.532, 4.544, 4.614, 4.617, 4.654, 4.680,
   4.743, 4.750, 4.787, 4.804, 4.839, 4.855, 4.877, 4.892,
   4.912, 4.962, 5.084, 5.113, 5.162, 5.060, 5.221, 5.202,
   5.251, 5.226, 5.312, 5.370, 5.342, 5.351, 5.367, 5.339,
   5.413, 5.402, 5.428, 5.436, 5.463, 5.476, 5.501, 5.521,
   5.526, 5.539, 5.655, 5.658, 5.684, 5.670, 5.710, 5.700,
   5.851, 5.744, 5.864, 5.905, 5.815, 5.815, 5.843, 5.850, 5.857};

static double _mserms[N_ELEMENTS] = {
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  3.0055    ,  2.9936    ,  3.057     ,  3.061     ,  3.1224    ,
  3.1889    ,  3.2611    ,  3.365     ,  3.4274    ,  3.4349    ,
  3.4776    ,  3.5459    ,  3.5921    ,  3.6002    ,  3.6452    ,
  3.7057    ,  3.7377    ,  3.7875    ,  3.7757    ,  3.8823    ,
  3.9283    ,  3.9973    ,  4.0742    ,  4.0968    ,  4.14      ,
  4.1629    ,  4.1884    ,  4.2036    ,  4.224     ,  4.243     ,
  4.2694    ,  4.324     ,  4.4091    ,  4.424     ,  4.4809    ,
  4.4945    ,  4.5318    ,  4.5454    ,  4.5944    ,  4.6156    ,
  4.6519    ,  4.6802    ,  4.7423    ,  4.75      ,  4.7859    ,
  4.8041    ,  4.8378    ,  4.855     ,  4.8771    ,  4.8919    ,
  4.9123    ,  4.91548125,  4.85956687,  4.88787188,  4.93128563,
  5.07265   ,  5.0898425 ,  5.17599   ,  5.238471  ,  5.303984  ,
  5.3108    ,  5.302875  ,  5.31358125,  5.31725813,  5.31381881,
  5.3230125 ,  5.331411  ,  5.37975   ,  5.41033488,  5.44389637,
  5.4648    ,  5.48103366,  5.49260437,  5.49895157,  5.51991148,
  5.53461255,  5.51322061,  5.54938223,  5.5893442 ,  5.65759688,
  5.6849025 ,  5.78639063,  5.8118025 ,  5.89141438,  5.88691415,
  5.90087899,  5.91494926,  5.91168879,  5.88264623,  5.86691016,
  5.85471211,  5.84942109,  5.86143996,  5.87414508,  5.89325189,
  5.919     ,  5.92795102,  5.94922305,  5.95825928,  5.97967266,
  5.99007373,  5.97927161,  5.97652411,  5.98987011,  6.02529002,
  6.08752438,  6.13926806,  6.15654237,  6.15927535,  6.16201645,
  6.17454773};
 
static double _errms[N_ELEMENTS][NISO];

static double _xfermi[NFERMI];
static double _yfermi[NFERMI];
static double _rfermi[5][NFERMI];

int InitNucleus() {
  _errms[0][0] = 0.8783;
  _errms[0][1] = 2.1421;
  _errms[0][2] = 1.7591;
  _errms[1][0] = 1.9661;
  _errms[1][1] = 1.6755;
  _errms[1][3] = 2.0660;
  _errms[1][5] = 1.9239;
  _errms[2][1] = 2.5890;
  _errms[2][2] = 2.4440;
  _errms[2][3] = 2.3390;
  _errms[2][4] = 2.2450;
  _errms[2][6] = 2.4820;
  _errms[3][0] = 2.6460;
  _errms[3][2] = 2.5190;
  _errms[3][3] = 2.3550;
  _errms[3][4] = 2.4630;
  _errms[4][1] = 2.4277;
  _errms[4][2] = 2.4060;
  _errms[5][1] = 2.4702;
  _errms[5][2] = 2.4614;
  _errms[5][3] = 2.5025;
  _errms[6][1] = 2.5582;
  _errms[6][2] = 2.6058;
  _errms[7][1] = 2.6991;
  _errms[7][2] = 2.6932;
  _errms[7][3] = 2.7726;
  _errms[8][2] = 2.8976;
  _errms[9][0] = 3.0082;
  _errms[9][1] = 3.0055;
  _errms[9][2] = 2.9695;
  _errms[9][3] = 2.9525;
  _errms[9][4] = 2.9104;
  _errms[9][5] = 2.9007;
  _errms[9][6] = 2.9316;
  _errms[9][7] = 2.9251;
  _errms[9][9] = 2.9642;
  _errms[9][56] = 3.0413;
  _errms[9][57] = 2.9714;
  _errms[10][0] = 3.0136;
  _errms[10][1] = 2.9852;
  _errms[10][2] = 2.9936;
  _errms[10][3] = 2.9735;
  _errms[10][4] = 2.9769;
  _errms[10][5] = 2.9928;
  _errms[10][6] = 3.0136;
  _errms[10][7] = 3.0400;
  _errms[10][8] = 3.0922;
  _errms[10][9] = 3.1180;
  _errms[10][10] = 3.1704;
  _errms[10][57] = 2.9718;
  _errms[11][1] = 3.0570;
  _errms[11][2] = 3.0284;
  _errms[11][3] = 3.0337;
  _errms[12][2] = 3.0610;
  _errms[13][1] = 3.1224;
  _errms[13][2] = 3.1176;
  _errms[13][3] = 3.1336;
  _errms[14][2] = 3.1889;
  _errms[15][1] = 3.2611;
  _errms[15][3] = 3.2847;
  _errms[15][5] = 3.2985;
  _errms[16][2] = 3.3654;
  _errms[16][4] = 3.3840;
  _errms[17][0] = 3.3636;
  _errms[17][1] = 3.3905;
  _errms[17][2] = 3.3908;
  _errms[17][3] = 3.4028;
  _errms[17][4] = 3.4093;
  _errms[17][5] = 3.4274;
  _errms[17][6] = 3.4251;
  _errms[17][7] = 3.4354;
  _errms[17][8] = 3.4414;
  _errms[17][9] = 3.4454;
  _errms[17][11] = 3.4377;
  _errms[17][55] = 3.3468;
  _errms[17][56] = 3.3438;
  _errms[17][57] = 3.3654;
  _errms[18][1] = 3.4264;
  _errms[18][2] = 3.4349;
  _errms[18][3] = 3.4381;
  _errms[18][4] = 3.4518;
  _errms[18][5] = 3.4517;
  _errms[18][6] = 3.4556;
  _errms[18][7] = 3.4563;
  _errms[18][8] = 3.4605;
  _errms[18][9] = 3.4558;
  _errms[18][10] = 3.4534;
  _errms[19][0] = 3.4595;
  _errms[19][1] = 3.4776;
  _errms[19][2] = 3.4780;
  _errms[19][3] = 3.5081;
  _errms[19][4] = 3.4954;
  _errms[19][5] = 3.5179;
  _errms[19][6] = 3.4944;
  _errms[19][7] = 3.4953;
  _errms[19][8] = 3.4783;
  _errms[19][9] = 3.4771;
  _errms[19][11] = 3.5168;
  _errms[20][1] = 3.5702;
  _errms[20][2] = 3.5575;
  _errms[20][3] = 3.5432;
  _errms[20][4] = 3.5459;
  _errms[20][5] = 3.5243;
  _errms[21][1] = 3.6115;
  _errms[21][2] = 3.5939;
  _errms[21][3] = 3.6070;
  _errms[21][4] = 3.5962;
  _errms[21][5] = 3.5921;
  _errms[21][6] = 3.5733;
  _errms[21][7] = 3.5704;
  _errms[22][6] = 3.6002;
  _errms[23][3] = 3.6588;
  _errms[23][5] = 3.6452;
  _errms[23][6] = 3.6511;
  _errms[23][7] = 3.6885;
  _errms[24][1] = 3.7120;
  _errms[24][2] = 3.7026;
  _errms[24][3] = 3.6706;
  _errms[24][4] = 3.6662;
  _errms[24][5] = 3.6834;
  _errms[24][6] = 3.7057;
  _errms[24][7] = 3.7146;
  _errms[25][3] = 3.6933;
  _errms[25][5] = 3.7377;
  _errms[25][6] = 3.7532;
  _errms[25][7] = 3.7745;
  _errms[26][6] = 3.7875;
  _errms[27][3] = 3.7757;
  _errms[27][5] = 3.8118;
  _errms[27][6] = 3.8225;
  _errms[27][7] = 3.8399;
  _errms[27][9] = 3.8572;
  _errms[28][6] = 3.8823;
  _errms[28][8] = 3.9022;
  _errms[29][5] = 3.9283;
  _errms[29][7] = 3.9491;
  _errms[29][8] = 3.9530;
  _errms[29][9] = 3.9658;
  _errms[29][11] = 3.9845;
  _errms[30][8] = 3.9973;
  _errms[30][10] = 4.0118;
  _errms[31][7] = 4.0414;
  _errms[31][9] = 4.0576;
  _errms[31][10] = 4.0632;
  _errms[31][11] = 4.0742;
  _errms[31][13] = 4.0811;
  _errms[32][10] = 4.0968;
  _errms[33][7] = 4.0700;
  _errms[33][9] = 4.1395;
  _errms[33][10] = 4.1395;
  _errms[33][11] = 4.1406;
  _errms[33][13] = 4.1400;
  _errms[33][15] = 4.1400;
  _errms[34][10] = 4.1629;
  _errms[34][12] = 4.1599;
  _errms[35][1] = 4.1635;
  _errms[35][3] = 4.1870;
  _errms[35][4] = 4.2097;
  _errms[35][5] = 4.2020;
  _errms[35][6] = 4.2082;
  _errms[35][7] = 4.2038;
  _errms[35][8] = 4.2034;
  _errms[35][9] = 4.1970;
  _errms[35][10] = 4.1952;
  _errms[35][11] = 4.1919;
  _errms[35][12] = 4.1871;
  _errms[35][13] = 4.1884;
  _errms[35][14] = 4.1846;
  _errms[35][15] = 4.1835;
  _errms[35][16] = 4.1984;
  _errms[35][17] = 4.2171;
  _errms[35][18] = 4.2286;
  _errms[35][19] = 4.2423;
  _errms[35][20] = 4.2543;
  _errms[35][21] = 4.2724;
  _errms[35][22] = 4.2794;
  _errms[35][23] = 4.3002;
  _errms[35][24] = 4.3067;
  _errms[35][25] = 4.3267;
  _errms[36][3] = 4.2273;
  _errms[36][4] = 4.2356;
  _errms[36][5] = 4.2385;
  _errms[36][6] = 4.2284;
  _errms[36][7] = 4.2271;
  _errms[36][8] = 4.2213;
  _errms[36][9] = 4.2160;
  _errms[36][10] = 4.2058;
  _errms[36][11] = 4.1999;
  _errms[36][12] = 4.2036;
  _errms[36][13] = 4.2025;
  _errms[36][14] = 4.1989;
  _errms[36][15] = 4.2170;
  _errms[36][16] = 4.2391;
  _errms[36][17] = 4.2554;
  _errms[36][18] = 4.2723;
  _errms[36][19] = 4.2903;
  _errms[36][20] = 4.3048;
  _errms[36][21] = 4.3184;
  _errms[36][22] = 4.3391;
  _errms[36][23] = 4.3501;
  _errms[36][24] = 4.4231;
  _errms[36][25] = 4.4336;
  _errms[37][2] = 4.2569;
  _errms[37][3] = 4.2561;
  _errms[37][4] = 4.2586;
  _errms[37][5] = 4.2562;
  _errms[37][6] = 4.2547;
  _errms[37][7] = 4.2478;
  _errms[37][8] = 4.2455;
  _errms[37][9] = 4.2394;
  _errms[37][10] = 4.2304;
  _errms[37][11] = 4.2307;
  _errms[37][12] = 4.2249;
  _errms[37][13] = 4.2240;
  _errms[37][14] = 4.2407;
  _errms[37][15] = 4.2611;
  _errms[37][16] = 4.2740;
  _errms[37][17] = 4.2924;
  _errms[37][18] = 4.3026;
  _errms[37][19] = 4.3191;
  _errms[37][20] = 4.3305;
  _errms[37][21] = 4.3522;
  _errms[37][22] = 4.3625;
  _errms[37][23] = 4.4377;
  _errms[37][24] = 4.4495;
  _errms[37][25] = 4.4640;
  _errms[38][9] = 4.2513;
  _errms[38][10] = 4.2498;
  _errms[38][11] = 4.2441;
  _errms[38][12] = 4.2430;
  _errms[38][13] = 4.2573;
  _errms[38][15] = 4.2887;
  _errms[38][16] = 4.3052;
  _errms[38][17] = 4.3142;
  _errms[38][18] = 4.3284;
  _errms[38][19] = 4.3402;
  _errms[38][20] = 4.3580;
  _errms[38][21] = 4.3711;
  _errms[38][22] = 4.4658;
  _errms[38][23] = 4.4705;
  _errms[38][24] = 4.4863;
  _errms[38][25] = 4.4911;
  _errms[39][8] = 4.2789;
  _errms[39][9] = 4.2787;
  _errms[39][10] = 4.2706;
  _errms[39][11] = 4.2694;
  _errms[39][12] = 4.2845;
  _errms[39][13] = 4.3057;
  _errms[39][15] = 4.3320;
  _errms[39][17] = 4.3512;
  _errms[39][18] = 4.3792;
  _errms[39][19] = 4.4012;
  _errms[39][20] = 4.4156;
  _errms[39][21] = 4.4891;
  _errms[39][22] = 4.5119;
  _errms[39][23] = 4.5292;
  _errms[40][9] = 4.2891;
  _errms[40][10] = 4.2878;
  _errms[40][11] = 4.3026;
  _errms[40][12] = 4.3240;
  _errms[40][18] = 4.4062;
  _errms[40][20] = 4.4861;
  _errms[40][22] = 4.5097;
  _errms[41][7] = 4.3265;
  _errms[41][8] = 4.3182;
  _errms[41][9] = 4.3151;
  _errms[41][11] = 4.3529;
  _errms[41][12] = 4.3628;
  _errms[41][13] = 4.3847;
  _errms[41][14] = 4.3880;
  _errms[41][15] = 4.4091;
  _errms[41][17] = 4.4468;
  _errms[41][19] = 4.4914;
  _errms[41][20] = 4.5145;
  _errms[41][21] = 4.5249;
  _errms[41][22] = 4.5389;
  _errms[41][23] = 4.5490;
  _errms[41][25] = 4.5602;
  _errms[43][9] = 4.3908;
  _errms[43][11] = 4.4229;
  _errms[43][12] = 4.4338;
  _errms[43][13] = 4.4531;
  _errms[43][14] = 4.4606;
  _errms[43][15] = 4.4809;
  _errms[43][17] = 4.5098;
  _errms[44][14] = 4.4945;
  _errms[45][11] = 4.4827;
  _errms[45][13] = 4.5078;
  _errms[45][14] = 4.5150;
  _errms[45][15] = 4.5318;
  _errms[45][17] = 4.5563;
  _errms[45][19] = 4.5782;
  _errms[46][8] = 4.4799;
  _errms[46][10] = 4.5036;
  _errms[46][11] = 4.5119;
  _errms[46][12] = 4.5269;
  _errms[46][14] = 4.5454;
  _errms[46][16] = 4.5638;
  _errms[47][7] = 4.4810;
  _errms[47][8] = 4.4951;
  _errms[47][9] = 4.5122;
  _errms[47][10] = 4.5216;
  _errms[47][11] = 4.5383;
  _errms[47][12] = 4.5466;
  _errms[47][13] = 4.5577;
  _errms[47][14] = 4.5601;
  _errms[47][15] = 4.5765;
  _errms[47][16] = 4.5845;
  _errms[47][17] = 4.5944;
  _errms[47][18] = 4.6012;
  _errms[47][19] = 4.6087;
  _errms[47][20] = 4.6114;
  _errms[47][21] = 4.6203;
  _errms[47][22] = 4.6136;
  _errms[47][23] = 4.6246;
  _errms[47][25] = 4.6300;
  _errms[48][7] = 4.5184;
  _errms[48][8] = 4.5311;
  _errms[48][9] = 4.5375;
  _errms[48][10] = 4.5494;
  _errms[48][11] = 4.5571;
  _errms[48][12] = 4.5685;
  _errms[48][13] = 4.5742;
  _errms[48][14] = 4.5856;
  _errms[48][15] = 4.5907;
  _errms[48][16] = 4.6010;
  _errms[48][17] = 4.6056;
  _errms[48][18] = 4.6156;
  _errms[48][19] = 4.6211;
  _errms[48][20] = 4.6292;
  _errms[48][21] = 4.6335;
  _errms[48][22] = 4.6407;
  _errms[48][23] = 4.6443;
  _errms[48][24] = 4.6505;
  _errms[48][25] = 4.6534;
  _errms[48][26] = 4.6594;
  _errms[48][27] = 4.6625;
  _errms[48][28] = 4.6670;
  _errms[48][29] = 4.6702;
  _errms[48][30] = 4.6733;
  _errms[49][9] = 4.5605;
  _errms[49][10] = 4.5679;
  _errms[49][11] = 4.5785;
  _errms[49][12] = 4.5836;
  _errms[49][13] = 4.5948;
  _errms[49][14] = 4.6015;
  _errms[49][15] = 4.6099;
  _errms[49][16] = 4.6148;
  _errms[49][17] = 4.6250;
  _errms[49][18] = 4.6302;
  _errms[49][19] = 4.6393;
  _errms[49][20] = 4.6438;
  _errms[49][21] = 4.6519;
  _errms[49][22] = 4.6566;
  _errms[49][23] = 4.6634;
  _errms[49][24] = 4.6665;
  _errms[49][25] = 4.6735;
  _errms[49][26] = 4.6765;
  _errms[49][27] = 4.6833;
  _errms[49][28] = 4.6867;
  _errms[49][29] = 4.6921;
  _errms[49][30] = 4.6934;
  _errms[49][31] = 4.7019;
  _errms[49][32] = 4.7078;
  _errms[49][33] = 4.7093;
  _errms[50][20] = 4.6802;
  _errms[50][22] = 4.6879;
  _errms[51][13] = 4.6847;
  _errms[51][15] = 4.6956;
  _errms[51][17] = 4.7038;
  _errms[51][19] = 4.7095;
  _errms[51][20] = 4.7117;
  _errms[51][21] = 4.7183;
  _errms[51][22] = 4.7204;
  _errms[51][23] = 4.7266;
  _errms[51][25] = 4.7346;
  _errms[51][27] = 4.7423;
  _errms[51][29] = 4.7500;
  _errms[51][31] = 4.7569;
  _errms[51][33] = 4.7815;
  _errms[52][22] = 4.7500;
  _errms[53][9] = 4.7211;
  _errms[53][11] = 4.7387;
  _errms[53][13] = 4.7509;
  _errms[53][15] = 4.7590;
  _errms[53][17] = 4.7661;
  _errms[53][19] = 4.7722;
  _errms[53][20] = 4.7747;
  _errms[53][21] = 4.7774;
  _errms[53][22] = 4.7775;
  _errms[53][23] = 4.7818;
  _errms[53][24] = 4.7808;
  _errms[53][25] = 4.7859;
  _errms[53][26] = 4.7831;
  _errms[53][27] = 4.7899;
  _errms[53][29] = 4.7964;
  _errms[53][30] = 4.8094;
  _errms[53][31] = 4.8279;
  _errms[53][32] = 4.8409;
  _errms[53][33] = 4.8566;
  _errms[53][34] = 4.8694;
  _errms[53][35] = 4.8841;
  _errms[53][36] = 4.8942;
  _errms[53][37] = 4.9082;
  _errms[53][39] = 4.9315;
  _errms[54][9] = 4.7832;
  _errms[54][10] = 4.7896;
  _errms[54][11] = 4.7915;
  _errms[54][12] = 4.7769;
  _errms[54][13] = 4.7773;
  _errms[54][14] = 4.7820;
  _errms[54][15] = 4.7828;
  _errms[54][16] = 4.7880;
  _errms[54][17] = 4.7872;
  _errms[54][18] = 4.7936;
  _errms[54][19] = 4.7921;
  _errms[54][20] = 4.7981;
  _errms[54][21] = 4.7992;
  _errms[54][22] = 4.8026;
  _errms[54][23] = 4.8002;
  _errms[54][24] = 4.8041;
  _errms[54][25] = 4.8031;
  _errms[54][26] = 4.8067;
  _errms[54][27] = 4.8059;
  _errms[54][28] = 4.8128;
  _errms[54][29] = 4.8255;
  _errms[54][30] = 4.8422;
  _errms[54][31] = 4.8554;
  _errms[54][32] = 4.8689;
  _errms[54][33] = 4.8825;
  _errms[54][34] = 4.8965;
  _errms[54][35] = 4.9055;
  _errms[54][36] = 4.9188;
  _errms[54][37] = 4.9281;
  _errms[55][9] = 4.8092;
  _errms[55][10] = 4.8176;
  _errms[55][11] = 4.8153;
  _errms[55][12] = 4.8135;
  _errms[55][13] = 4.8185;
  _errms[55][14] = 4.8177;
  _errms[55][15] = 4.8221;
  _errms[55][16] = 4.8204;
  _errms[55][17] = 4.8255;
  _errms[55][18] = 4.8248;
  _errms[55][19] = 4.8283;
  _errms[55][20] = 4.8276;
  _errms[55][21] = 4.8303;
  _errms[55][22] = 4.8286;
  _errms[55][23] = 4.8322;
  _errms[55][24] = 4.8294;
  _errms[55][25] = 4.8334;
  _errms[55][26] = 4.8314;
  _errms[55][27] = 4.8378;
  _errms[55][28] = 4.8513;
  _errms[55][29] = 4.8684;
  _errms[55][30] = 4.8807;
  _errms[55][31] = 4.8953;
  _errms[55][32] = 4.9087;
  _errms[55][33] = 4.9236;
  _errms[55][34] = 4.9345;
  _errms[55][35] = 4.9479;
  _errms[55][37] = 4.9731;
  _errms[56][22] = 4.8488;
  _errms[56][24] = 4.8496;
  _errms[56][25] = 4.8473;
  _errms[56][26] = 4.8550;
  _errms[57][21] = 4.8739;
  _errms[57][23] = 4.8737;
  _errms[57][25] = 4.8771;
  _errms[57][27] = 4.9063;
  _errms[57][29] = 4.9303;
  _errms[57][31] = 4.9590;
  _errms[57][33] = 4.9893;
  _errms[58][24] = 4.8919;
  _errms[59][13] = 4.9174;
  _errms[59][15] = 4.9128;
  _errms[59][16] = 4.9086;
  _errms[59][17] = 4.9111;
  _errms[59][18] = 4.9080;
  _errms[59][19] = 4.9123;
  _errms[59][20] = 4.9076;
  _errms[59][21] = 4.9101;
  _errms[59][22] = 4.9057;
  _errms[59][23] = 4.9123;
  _errms[59][24] = 4.9254;
  _errms[59][25] = 4.9421;
  _errms[59][26] = 4.9535;
  _errms[59][27] = 4.9696;
  _errms[59][29] = 4.9999;
  _errms[59][31] = 5.0400;
  _errms[61][15] = 4.9599;
  _errms[61][16] = 4.9556;
  _errms[61][17] = 4.9565;
  _errms[61][18] = 4.9517;
  _errms[61][19] = 4.9518;
  _errms[61][20] = 4.9479;
  _errms[61][21] = 4.9524;
  _errms[61][22] = 4.9651;
  _errms[61][23] = 4.9808;
  _errms[61][24] = 4.9892;
  _errms[61][25] = 5.0042;
  _errms[61][26] = 5.0134;
  _errms[61][27] = 5.0387;
  _errms[61][28] = 5.0550;
  _errms[61][29] = 5.0819;
  _errms[61][30] = 5.0925;
  _errms[61][31] = 5.1053;
  _errms[62][12] = 4.9762;
  _errms[62][13] = 4.9779;
  _errms[62][14] = 4.9760;
  _errms[62][15] = 4.9695;
  _errms[62][16] = 4.9697;
  _errms[62][17] = 4.9607;
  _errms[62][18] = 4.9636;
  _errms[62][19] = 4.9612;
  _errms[62][20] = 4.9663;
  _errms[62][21] = 4.9789;
  _errms[62][22] = 4.9938;
  _errms[62][23] = 5.0045;
  _errms[62][24] = 5.0202;
  _errms[62][25] = 5.0296;
  _errms[62][26] = 5.0522;
  _errms[62][27] = 5.1064;
  _errms[62][28] = 5.1115;
  _errms[62][29] = 5.1239;
  _errms[62][30] = 5.1221;
  _errms[62][31] = 5.1264;
  _errms[62][32] = 5.1351;
  _errms[62][33] = 5.1413;
  _errms[62][34] = 5.1498;
  _errms[63][18] = 4.9786;
  _errms[63][19] = 4.9801;
  _errms[63][21] = 5.0080;
  _errms[63][23] = 5.0342;
  _errms[63][25] = 5.0774;
  _errms[63][27] = 5.1223;
  _errms[63][28] = 5.1319;
  _errms[63][29] = 5.1420;
  _errms[63][30] = 5.1449;
  _errms[63][31] = 5.1569;
  _errms[63][33] = 5.1734;
  _errms[64][18] = 4.9201;
  _errms[64][19] = 4.9291;
  _errms[64][20] = 4.9427;
  _errms[64][21] = 4.9499;
  _errms[64][22] = 4.9630;
  _errms[64][23] = 4.9689;
  _errms[64][24] = 4.9950;
  _errms[64][25] = 5.0333;
  _errms[64][26] = 5.0391;
  _errms[64][28] = 5.0489;
  _errms[64][30] = 5.0600;
  _errms[65][15] = 5.0438;
  _errms[65][17] = 5.0455;
  _errms[65][18] = 5.0567;
  _errms[65][19] = 5.0706;
  _errms[65][20] = 5.0801;
  _errms[65][21] = 5.0950;
  _errms[65][22] = 5.1035;
  _errms[65][23] = 5.1241;
  _errms[65][24] = 5.1457;
  _errms[65][25] = 5.1622;
  _errms[65][26] = 5.1709;
  _errms[65][27] = 5.1815;
  _errms[65][28] = 5.1825;
  _errms[65][29] = 5.1951;
  _errms[65][30] = 5.1962;
  _errms[65][31] = 5.2074;
  _errms[65][32] = 5.2099;
  _errms[65][33] = 5.2218;
  _errms[66][18] = 5.0398;
  _errms[66][19] = 5.0614;
  _errms[66][20] = 5.0760;
  _errms[66][21] = 5.0856;
  _errms[66][22] = 5.1076;
  _errms[66][23] = 5.1156;
  _errms[66][24] = 5.1535;
  _errms[66][25] = 5.1571;
  _errms[66][26] = 5.1675;
  _errms[66][27] = 5.1662;
  _errms[66][28] = 5.1785;
  _errms[66][29] = 5.1817;
  _errms[66][30] = 5.1907;
  _errms[66][32] = 5.2022;
  _errms[67][15] = 5.0548;
  _errms[67][17] = 5.0843;
  _errms[67][19] = 5.1129;
  _errms[67][21] = 5.1429;
  _errms[67][23] = 5.1761;
  _errms[67][25] = 5.2045;
  _errms[67][27] = 5.2246;
  _errms[67][29] = 5.2389;
  _errms[67][31] = 5.2516;
  _errms[67][32] = 5.2560;
  _errms[67][33] = 5.2644;
  _errms[67][35] = 5.2789;
  _errms[68][16] = 5.0643;
  _errms[68][17] = 5.0755;
  _errms[68][19] = 5.0976;
  _errms[68][20] = 5.1140;
  _errms[68][21] = 5.1235;
  _errms[68][22] = 5.1392;
  _errms[68][23] = 5.1504;
  _errms[68][24] = 5.1616;
  _errms[68][25] = 5.1713;
  _errms[68][26] = 5.1849;
  _errms[68][27] = 5.1906;
  _errms[68][28] = 5.2004;
  _errms[68][29] = 5.2046;
  _errms[68][30] = 5.2129;
  _errms[68][31] = 5.2170;
  _errms[68][32] = 5.2256;
  _errms[68][33] = 5.2303;
  _errms[68][34] = 5.2388;
  _errms[68][35] = 5.2411;
  _errms[69][13] = 5.0423;
  _errms[69][15] = 5.0875;
  _errms[69][16] = 5.1040;
  _errms[69][17] = 5.1219;
  _errms[69][18] = 5.1324;
  _errms[69][19] = 5.1498;
  _errms[69][20] = 5.1629;
  _errms[69][21] = 5.1781;
  _errms[69][22] = 5.1889;
  _errms[69][23] = 5.2054;
  _errms[69][24] = 5.2157;
  _errms[69][25] = 5.2307;
  _errms[69][26] = 5.2399;
  _errms[69][27] = 5.2525;
  _errms[69][28] = 5.2621;
  _errms[69][29] = 5.2702;
  _errms[69][30] = 5.2771;
  _errms[69][31] = 5.2853;
  _errms[69][32] = 5.2906;
  _errms[69][33] = 5.2995;
  _errms[69][34] = 5.3046;
  _errms[69][35] = 5.3108;
  _errms[69][36] = 5.3135;
  _errms[69][37] = 5.3215;
  _errms[70][20] = 5.2293;
  _errms[70][21] = 5.2398;
  _errms[70][22] = 5.2567;
  _errms[70][23] = 5.2677;
  _errms[70][24] = 5.2830;
  _errms[70][25] = 5.2972;
  _errms[70][26] = 5.3108;
  _errms[70][27] = 5.3227;
  _errms[70][28] = 5.3290;
  _errms[70][29] = 5.3364;
  _errms[70][30] = 5.3436;
  _errms[70][31] = 5.3486;
  _errms[70][32] = 5.3577;
  _errms[70][33] = 5.3634;
  _errms[70][34] = 5.3700;
  _errms[70][35] = 5.3739;
  _errms[70][36] = 5.3815;
  _errms[70][37] = 5.3857;
  _errms[70][38] = 5.3917;
  _errms[71][27] = 5.2898;
  _errms[71][28] = 5.3041;
  _errms[71][29] = 5.3065;
  _errms[71][30] = 5.3140;
  _errms[71][31] = 5.3201;
  _errms[71][32] = 5.3191;
  _errms[71][33] = 5.3286;
  _errms[71][34] = 5.3309;
  _errms[71][35] = 5.3371;
  _errms[71][36] = 5.3408;
  _errms[71][37] = 5.3470;
  _errms[71][39] = 5.3516;
  _errms[72][36] = 5.3507;
  _errms[73][33] = 5.3491;
  _errms[73][35] = 5.3559;
  _errms[73][36] = 5.3611;
  _errms[73][37] = 5.3658;
  _errms[73][39] = 5.3743;
  _errms[74][36] = 5.3596;
  _errms[74][38] = 5.3698;
  _errms[75][33] = 5.3823;
  _errms[75][35] = 5.3909;
  _errms[75][36] = 5.3933;
  _errms[75][37] = 5.3993;
  _errms[75][38] = 5.4016;
  _errms[75][39] = 5.4062;
  _errms[75][41] = 5.4126;
  _errms[76][29] = 5.3705;
  _errms[76][30] = 5.3780;
  _errms[76][31] = 5.3805;
  _errms[76][32] = 5.3854;
  _errms[76][33] = 5.3900;
  _errms[76][34] = 5.3812;
  _errms[76][35] = 5.3838;
  _errms[76][36] = 5.3898;
  _errms[76][38] = 5.3968;
  _errms[76][40] = 5.4032;
  _errms[77][23] = 5.3728;
  _errms[77][24] = 5.3915;
  _errms[77][25] = 5.3891;
  _errms[77][26] = 5.3996;
  _errms[77][27] = 5.3969;
  _errms[77][28] = 5.4038;
  _errms[77][29] = 5.4015;
  _errms[77][30] = 5.4148;
  _errms[77][31] = 5.4037;
  _errms[77][32] = 5.4063;
  _errms[77][33] = 5.4053;
  _errms[77][34] = 5.4060;
  _errms[77][35] = 5.4108;
  _errms[77][36] = 5.4102;
  _errms[77][37] = 5.4169;
  _errms[77][38] = 5.4191;
  _errms[77][39] = 5.4236;
  _errms[77][40] = 5.4270;
  _errms[77][41] = 5.4307;
  _errms[77][43] = 5.4383;
  _errms[78][26] = 5.4247;
  _errms[78][27] = 5.4306;
  _errms[78][28] = 5.4296;
  _errms[78][29] = 5.4354;
  _errms[78][30] = 5.4018;
  _errms[78][31] = 5.4049;
  _errms[78][32] = 5.4084;
  _errms[78][33] = 5.4109;
  _errms[78][34] = 5.4147;
  _errms[78][35] = 5.4179;
  _errms[78][36] = 5.4221;
  _errms[78][37] = 5.4252;
  _errms[78][38] = 5.4298;
  _errms[78][39] = 5.4332;
  _errms[78][40] = 5.4371;
  _errms[78][41] = 5.4400;
  _errms[78][42] = 5.4454;
  _errms[79][22] = 5.4364;
  _errms[79][23] = 5.3833;
  _errms[79][24] = 5.4405;
  _errms[79][25] = 5.3949;
  _errms[79][26] = 5.4397;
  _errms[79][27] = 5.4017;
  _errms[79][28] = 5.4046;
  _errms[79][29] = 5.4085;
  _errms[79][30] = 5.4100;
  _errms[79][31] = 5.4158;
  _errms[79][32] = 5.4171;
  _errms[79][33] = 5.4232;
  _errms[79][34] = 5.4238;
  _errms[79][35] = 5.4309;
  _errms[79][36] = 5.4345;
  _errms[79][37] = 5.4385;
  _errms[79][38] = 5.4412;
  _errms[79][39] = 5.4463;
  _errms[79][40] = 5.4474;
  _errms[79][41] = 5.4551;
  _errms[79][42] = 5.4581;
  _errms[79][43] = 5.4648;
  _errms[79][44] = 5.4679;
  _errms[79][45] = 5.4744;
  _errms[79][46] = 5.4776;
  _errms[79][47] = 5.4837;
  _errms[80][27] = 5.4017;
  _errms[80][29] = 5.4121;
  _errms[80][30] = 5.4169;
  _errms[80][31] = 5.4191;
  _errms[80][32] = 5.4243;
  _errms[80][33] = 5.4259;
  _errms[80][34] = 5.4325;
  _errms[80][35] = 5.4327;
  _errms[80][36] = 5.4388;
  _errms[80][37] = 5.4396;
  _errms[80][38] = 5.4479;
  _errms[80][39] = 5.4491;
  _errms[80][40] = 5.4573;
  _errms[80][41] = 5.4595;
  _errms[80][42] = 5.4666;
  _errms[80][43] = 5.4704;
  _errms[80][44] = 5.4759;
  _errms[80][46] = 5.4853;
  _errms[80][47] = 5.4946;
  _errms[81][19] = 5.3788;
  _errms[81][20] = 5.3869;
  _errms[81][21] = 5.3930;
  _errms[81][22] = 5.3984;
  _errms[81][23] = 5.4027;
  _errms[81][24] = 5.4079;
  _errms[81][25] = 5.4139;
  _errms[81][26] = 5.4177;
  _errms[81][27] = 5.4222;
  _errms[81][28] = 5.4229;
  _errms[81][29] = 5.4300;
  _errms[81][30] = 5.4310;
  _errms[81][31] = 5.4372;
  _errms[81][32] = 5.4389;
  _errms[81][33] = 5.4444;
  _errms[81][34] = 5.4446;
  _errms[81][35] = 5.4524;
  _errms[81][36] = 5.4529;
  _errms[81][37] = 5.4611;
  _errms[81][38] = 5.4629;
  _errms[81][39] = 5.4705;
  _errms[81][40] = 5.4727;
  _errms[81][41] = 5.4803;
  _errms[81][42] = 5.4828;
  _errms[81][43] = 5.4902;
  _errms[81][44] = 5.4943;
  _errms[81][45] = 5.5012;
  _errms[81][46] = 5.5100;
  _errms[81][47] = 5.5208;
  _errms[81][48] = 5.5290;
  _errms[81][49] = 5.5396;
  _errms[81][51] = 5.5577;
  _errms[82][37] = 5.4840;
  _errms[82][38] = 5.4911;
  _errms[82][39] = 5.4934;
  _errms[82][40] = 5.5008;
  _errms[82][41] = 5.5034;
  _errms[82][42] = 5.5103;
  _errms[82][43] = 5.5147;
  _errms[82][44] = 5.5211;
  _errms[82][45] = 5.5300;
  _errms[82][47] = 5.5489;
  _errms[82][48] = 5.5586;
  _errms[83][25] = 5.5220;
  _errms[83][27] = 5.5167;
  _errms[83][29] = 5.5136;
  _errms[83][31] = 5.5146;
  _errms[83][33] = 5.5199;
  _errms[83][35] = 5.5281;
  _errms[83][37] = 5.5378;
  _errms[83][38] = 5.5389;
  _errms[83][39] = 5.5480;
  _errms[83][40] = 5.5501;
  _errms[83][41] = 5.5584;
  _errms[83][42] = 5.5628;
  _errms[83][43] = 5.5704;
  _errms[83][49] = 5.6359;
  _errms[83][51] = 5.6558;
  _errms[85][31] = 5.5521;
  _errms[85][33] = 5.5568;
  _errms[85][34] = 5.5569;
  _errms[85][35] = 5.5640;
  _errms[85][36] = 5.5652;
  _errms[85][37] = 5.5725;
  _errms[85][38] = 5.5743;
  _errms[85][39] = 5.5813;
  _errms[85][40] = 5.5850;
  _errms[85][41] = 5.5915;
  _errms[85][47] = 5.6540;
  _errms[85][48] = 5.6648;
  _errms[85][49] = 5.6731;
  _errms[85][50] = 5.6834;
  _errms[85][51] = 5.6915;
  _errms[86][34] = 5.5720;
  _errms[86][35] = 5.5729;
  _errms[86][36] = 5.5799;
  _errms[86][37] = 5.5818;
  _errms[86][38] = 5.5882;
  _errms[86][39] = 5.5915;
  _errms[86][40] = 5.5977;
  _errms[86][47] = 5.6688;
  _errms[86][48] = 5.6790;
  _errms[86][49] = 5.6890;
  _errms[86][50] = 5.6951;
  _errms[86][51] = 5.7061;
  _errms[86][52] = 5.7112;
  _errms[86][53] = 5.7190;
  _errms[86][54] = 5.7335;
  _errms[86][55] = 5.7399;
  _errms[87][33] = 5.5850;
  _errms[87][34] = 5.5853;
  _errms[87][35] = 5.5917;
  _errms[87][36] = 5.5929;
  _errms[87][37] = 5.5991;
  _errms[87][38] = 5.6020;
  _errms[87][39] = 5.6079;
  _errms[87][45] = 5.6683;
  _errms[87][46] = 5.6795;
  _errms[87][47] = 5.6874;
  _errms[87][48] = 5.6973;
  _errms[87][49] = 5.7046;
  _errms[87][50] = 5.7150;
  _errms[87][51] = 5.7211;
  _errms[87][52] = 5.7283;
  _errms[87][53] = 5.7370;
  _errms[87][54] = 5.7455;
  _errms[87][55] = 5.7551;
  _errms[87][57] = 5.7714;
  _errms[89][48] = 5.7404;
  _errms[89][49] = 5.7488;
  _errms[89][50] = 5.7557;
  _errms[89][51] = 5.7670;
  _errms[89][53] = 5.7848;
  _errms[91][50] = 5.8203;
  _errms[91][51] = 5.8291;
  _errms[91][52] = 5.8337;
  _errms[91][53] = 5.8431;
  _errms[91][55] = 5.8571;
  _errms[93][51] = 5.8535;
  _errms[93][52] = 5.8601;
  _errms[93][53] = 5.8701;
  _errms[93][54] = 5.8748;
  _errms[93][55] = 5.8823;
  _errms[93][57] = 5.8948;
  _errms[94][52] = 5.8928;
  _errms[94][54] = 5.9048;
  _errms[95][51] = 5.8285;
  _errms[95][53] = 5.8429;
  _errms[95][54] = 5.8475;
  _errms[95][55] = 5.8562;
  _errms[95][57] = 5.8687;

  _xfermi[0] = -10.0;
  double dx = 40.0/(NFERMI-1);
  int i, k;
  for (i = 1; i < NFERMI; i++) {
    _xfermi[i] = _xfermi[i-1] + dx;
  }
  for (i = 0; i < NFERMI; i++) {
    _yfermi[i] = dx/(1 + exp(_xfermi[i]));
  }
  for (k = 0; k <= 4; k++) {
    _rfermi[k][0] = 0.0;
    NewtonCotes(_rfermi[k], _yfermi, 0, NFERMI-1, -1, 0);
    if (k < 4) {
      for (i = 0; i < NFERMI; i++) {
	_yfermi[i] *= _xfermi[i];
      }
    }
  }

  SetExtraPotential(-1, 0, NULL);
  return 0;
}

void SetExtraPotential(int m, int n, double *p) {
  int i;
  if (m < 0) {
    atom.nep = 0;
    for (i = 0; i < NEP; i++) {
      atom.epm[i] = m;      
    }
  } else {
    i = atom.nep;
    if (i == NEP) {
      printf("extra potential terms exceeded max: %d\n", NEP);
      return;
    }
    atom.epm[i] = m;
    if (m >= 0 && n > 0 && p != NULL) {
      memcpy(atom.epp[i], p, sizeof(double)*Min(n, NEPP));
    }
    atom.nep++;
  }
}

char *GetAtomicSymbolTable(void) {
  return (char *) _ename;
}

double *GetAtomicMassTable(void) {
  return _emass;
}

void IntegrateFermi(int nk, double *r, double x) {
  int k, one=1, np = 1, n = NFERMI;
  
  if (x >= _xfermi[0]) {
    for (k = 0; k < nk; k++) {
      UVIP3P(np, n, _xfermi, _rfermi[k], one, &x, &r[k]);
    }
  } else {
    double x1 = x;
    double x0 = _xfermi[0];
    r[0] = x1 - x0;
    if (nk > 1) {
      x1 *= x;
      x0 *= _xfermi[0];
      r[1] = (x1 - x0)/2.0;
      if (nk > 2) {
	x1 *= x;
	x0 *= _xfermi[0];
	r[2] = (x1 - x0)/3.0;
	if (nk > 3) {
	  x1 *= x;
	  x0 *= _xfermi[0];
	  r[3] = (x1 - x0)/4.0;
	  if (nk > 4) {
	    x1 *= x;
	    x0 *= _xfermi[0];
	    r[4] = (x1 - x0)/5.0;
	  }
	}
      }
    }
  }
}

double DiffRRMS(double c, double a, double a2, double a3, double a4, double r2,
		double y0, double y1, double y2, double y3, double y4) {
  double f, c2, c3, c4, r0, r1;
  c2 = c*c;
  c3 = c2*c;
  c4 = c3*c;
  IntegrateFermi(5, atom.rfermi, -c/a);
  y0 -= atom.rfermi[0];
  y1 -= atom.rfermi[1];
  y2 -= atom.rfermi[2];
  y3 -= atom.rfermi[3];
  y4 -= atom.rfermi[4];
  r0 = a2*y2 + 2*a*c*y1 + c2*y0;
  r1 = a4*y4 + 4*a3*c*y3 + 6*a2*c2*y2 + 4*a*c3*y1 + c4*y0;
  
  return r1/r0 - r2;
}

double GraspRRMS(double z, double m) {
  int iz, ia, i, k;
  double r0;
  iz = (int)(z-1);
  ia = (int)(1.5 + m - 2*(iz+1));
  if (ia >= 0 && ia < NISO) {
    r0 = _errms[iz][ia];
    if (r0 > 0) return r0;
  }
  for (i = 1; i < 10; i++) {
    k = ia-i;
    if (k >= 0 && k < NISO) {
      r0 = _errms[iz][k];
      if (r0 > 0) {
	r0 += 0.836*(pow(m,0.333)-pow(m-i,0.333));
	return r0;
      }
    }
    k = ia+i;
    if (k >= 0 && k < NISO) {
      r0 = _errms[iz][k];
      if (r0 > 0) {
	r0 += 0.836*(pow(m,0.333)-pow(m+i,0.333));
	return r0;
      }
    }
  }
  return 0.0;
}

int SetAtom(char *s, double z, double mass, double rn, double a, double rmse) {
  int i;
  char un[3] = "Xx";
  if (s == NULL || strlen(s) == 0) {
    if (z <= 0) {
      printf("atomic symbol and z cannot be both unset\n");
      return -1;
    }
    int iz = (int)z;
    if (iz <= N_ELEMENTS) {
      s = _ename[iz-1];
    } else {
      s = un;
    }
  } else {
    int iz = atoi(s);
    if (iz > 0) {      
      z = (double)iz;
      if (iz <= N_ELEMENTS) {
	s = _ename[iz-1];
      } else {
	s = un;
      }
    }
  }
  strncpy(atom.symbol, s, 2);
  atom.z0 = 0;
  atom.m0 = 0;
  for (i = 0; i < N_ELEMENTS; i++) {
    if (strncasecmp(_ename[i], s, 2) == 0) {
      atom.z0 = i+1;
      if (z <= 0) atom.atomic_number = i+1;
      atom.m0 = _emass[i];
      if (mass <= 0) atom.mass = _emass[i];
      break;
    }
  }
  if (z <= 0 || mass <= 0) {
    if (i == N_ELEMENTS) {
      printf("unknown element must have z and mass set\n");
      return -1;
    }
  }
  if (z > 0) {
    atom.atomic_number = z;
  } 
  if (mass > 0.0) {
    atom.mass = mass;
  }
  atom.rms0 = 1e-5*(0.570 + 0.836*pow(atom.mass,0.3333333333))/RBOHR;
  if (rn < 0.0) {
    if (rn > -1.5) {
      atom.rn = NucleusRRMS(atom.atomic_number);
    } else if (rn > -3.5) {
      if (atom.m0 > 0) {
	atom.rn = GraspRRMS(atom.atomic_number, atom.m0);
      } else {
	atom.rn = GraspRRMS(atom.atomic_number, atom.mass);
      }
      if (atom.rn <= 0) {
	if (rn > -2.5) {
	  atom.rn = NucleusRRMS(atom.atomic_number);
	} else {
	  atom.rn = 0.570 + 0.836*pow(atom.mass,0.3333333333);
	}
      }
    } else {
      int iz = (int)(atom.atomic_number-0.5);
      if (iz >= 0 && iz < N_ELEMENTS) {
	atom.rn = _arrms[iz];
      }
      if (atom.rn <= 0) {
	atom.rn = NucleusRRMS(atom.atomic_number);
      }
    }
    if (atom.rn <= 0) {
      atom.rn = 0.570 + 0.836*pow(atom.mass,0.3333333333);
    }
    if (atom.mass != atom.m0 && atom.m0 > 0) {
      atom.rn += 0.836*(pow(atom.mass,0.3333333)-pow(atom.m0,0.3333333));
    }
    atom.rn *= 1e-5/RBOHR;
    atom.rms = atom.rn;
  } else {
    rn *= 1e-5/RBOHR;
    atom.rn = rn;
    atom.rms = rn;
  }
  if (a >= 0) {
    atom.a = a;
  } else {
    atom.a = 2.3*1e-5/(RBOHR*4*log(3.0));
  }

  if (atom.a <= 0 && atom.rn > 0) {
    atom.rn = sqrt(5.0/3.0)*atom.rn;
    atom.z1 = 1.5*atom.atomic_number/atom.rn;
  }
  
  if (atom.rn > 0 && atom.a > 0) {
    i = NFERMI-1;
    a = atom.a;
    double y0 = _rfermi[0][i];
    double y1 = _rfermi[1][i];
    double y2 = _rfermi[2][i];
    double y3 = _rfermi[3][i];
    double y4 = _rfermi[4][i];
    double c0, c1, c, f;
    double a2 = a*a;
    double a3 = a2*a;
    double a4 = a2*a2;
    double r2 = atom.rn*atom.rn;
    c0 = 1e-2*atom.rn;
    c1 = 10*(atom.rn + atom.a);
    for (i = 0; i < 500; i++) {
      if (fabs(c1/c0-1) < 1e-8) break;      
      c = 0.5*(c0 + c1);
      f = DiffRRMS(c, a, a2, a3, a4, r2, y0, y1, y2, y3, y4);
      if (f < 0) {
	c0 = c;
      } else if (f > 0) {
	c1 = c;
      }	else {
	break;
      }
    }
    atom.c = c;
    if (i == 500) {
      printf("max iteration reached in determining fermi c param: %g %g %g %g %g %g %g %g %g %g %g\n",
	     atom.rms*1e5*RBOHR, atom.rn*1e5*RBOHR, atom.a, c1, c0, fabs(c1/c0-1), y0, y1, y2, y3, y4);
    }

    IntegrateFermi(3, atom.rfermi, -atom.c/atom.a);    
    atom.b = c*c*(y0-atom.rfermi[0]);
    atom.b += 2*a*c*(y1-atom.rfermi[1]);
    atom.b += a2*(y2-atom.rfermi[2]);  
    atom.b = atom.atomic_number/atom.b;
    atom.z1 = atom.c*(_rfermi[0][NFERMI-1]-atom.rfermi[0]);
    atom.z1 += atom.a*(_rfermi[1][NFERMI-1] - atom.rfermi[1]);
    atom.z1 *= atom.b;
  }

  if (rmse < 0) {
    int iz = (int)atom.atomic_number;
    if (rmse < 0 && iz <= N_ELEMENTS) {
      rmse = _mserms[iz-1];    
      if (rmse <= 0) {
	atom.rmse = atom.rms*1e5*RBOHR;
      } else {
	atom.rmse = rmse;
      }
    }
  } else {
    atom.rmse = rmse;
  }

  if (atom.atomic_number >= 10 && atom.atomic_number <= 120) {
    INIQED(atom.atomic_number, 9, atom.rn>0, atom.rmse);
  }
  return 0;
}

NUCLEUS *GetAtomicNucleus() {
  return &atom;
}

void PrintNucleus(int m, char *fn) {
  FILE *f;

  if (fn == NULL || strlen(fn) == 0 || strcmp(fn, "-")==0) {
    f = stdout;
  } else {
    f = fopen(fn, "w");
    if (f == NULL) {
      printf("cannot open iso output file: %s\n", fn);
      return;
    }
  }
  if (m == 0) {
    fprintf(f, "atom: %s\n", atom.symbol);
    fprintf(f, "z: %g\n", atom.atomic_number);
    fprintf(f, "mass: %g\n", atom.mass);
    fprintf(f, "rn: %g\n", atom.rn);
    fprintf(f, "rms: %g\n", atom.rms*1e5*RBOHR);
    fprintf(f, "rmse: %g\n", atom.rmse);
    fprintf(f, "fz1: %g\n", atom.z1);
    fprintf(f, "fa: %g\n", atom.a);
    fprintf(f, "fb: %g\n", atom.b);
    fprintf(f, "fc: %g\n", atom.c);
    int i;
    for (i = 0; i < atom.nep; i++) {
      fprintf(f, "ep: %d %d %g %g\n",
	      i, atom.epm[i], atom.epp[i][0], atom.epp[i][1]);
    }
  } else {
    fprintf(f, "Atomic number:\n");
    fprintf(f, "%15.8E\n", atom.atomic_number);
    fprintf(f, "Mass number (integer) :\n");
    fprintf(f, "%15.8E\n", (double)(int)(0.5+atom.mass));
    fprintf(f, "Fermi distribution parameter a:\n");
    fprintf(f, "%15.8E\n", atom.a*1e5*RBOHR);
    fprintf(f, "Fermi distribution parameter c:\n");
    fprintf(f, "%15.8E\n", atom.c*1e5*RBOHR);
    fprintf(f, "Mass of nucleus (in amu):\n");
    fprintf(f, "%15.8E\n", atom.mass);
    fprintf(f, "Nuclear spin (I) (in units of h / 2 pi):\n");
    fprintf(f, "%15.8E\n", 0.0);
    fprintf(f, "Nuclear dipole moment (in nuclear magnetons):\n");
    fprintf(f, "%15.8E\n", 0.0);
    fprintf(f, "Nuclear quadrupole moment (in barns):\n");
    fprintf(f, "%15.8E\n", 0.0);
    fprintf(f, "RNT:\n");
    fprintf(f, "%15.8E\n", 0.0);
    fprintf(f, "H:\n");
    fprintf(f, "%15.8E\n", 0.0);
    fprintf(f, "HP:\n");
    fprintf(f, "%15.8E\n", 0.0);
    fprintf(f, "NNNP:\n");
    fprintf(f, "%d\n", 0);
  }
  fflush(f);
  if (f != stdout) {
    fclose(f);
  }
}

double GetAtomicMass(void) {
  return atom.mass;
}

double GetAtomicNumber(void) {
  return atom.atomic_number;
}

char *GetAtomicSymbol(void) {
  return atom.symbol;
}

double GetAtomicChargeDist(double r) {
  double x;
  double r3;
  
  if (atom.rn <= 0) return 0.0;
  if (atom.a <= 0) {
    if (r > atom.rn) return 0.0;
    r3 = atom.rn;
    r3 = r3*r3*r3;
    return 3*atom.atomic_number/r3;
  }

  x = (r-atom.c)/atom.a;
  if (x < -10) r3 = atom.b;
  else if (x > 10) r3 = atom.b*exp(-x);
  else r3 = atom.b/(1+exp(x));
  r3 /= atom.a;
  return r3;
}

double GetExtraZ(double r, int i) {
  double z, zi;
  if (i >= atom.nep) return 0;
  z = 0.0;
  switch (atom.epm[i]) {
  case 0:
  case 100:    
    zi = atom.epp[i][0]*(atom.mass-atom.atomic_number);
    if (atom.epp[i][1] > 0) {
      zi *= exp(-2.68e-4*atom.epp[i][1]*r);
    }
    z += zi;
    break;
  case 1:
  case 101:
    if (r < atom.rms*atom.epp[i][1]) r = atom.rms*atom.epp[i][1];
    zi = 1.48e-8*(0.76+2.79/pow(atom.mass,0.33333))*atom.mass;
    zi *= atom.epp[i][0]*atom.rms*atom.rms/(r*r*r);
    zi *= FINE_STRUCTURE_CONST;
    z += zi;
    break;
  default:
    break;
  }
  return z;
}

double GetAtomicEffectiveZ(double r) {
  double x, y[3], z;
  int np = 3;
  int n = NFERMI;
  int one = 1;

  if (atom.rn <= 0) return (double)(atom.atomic_number);
  if (atom.a <= 0) {    
    if (r > atom.rn) {
      return (double) atom.atomic_number;
    } else {
      x = r/atom.rn;
      z = 3.0 - x*x;
      z = x*z*0.5*(atom.atomic_number);
      return z;
    }
  }
  x = (r - atom.c)/atom.a;
  if (x >= _xfermi[n-1]) {
    return (double)(atom.atomic_number);
  } else {
    IntegrateFermi(3, y, x);
    z = atom.c*atom.c*(y[0] - atom.rfermi[0]);
    z += 2*atom.c*atom.a*(y[1] - atom.rfermi[1]);
    z += atom.a*atom.a*(y[2] - atom.rfermi[2]);
    z += r*(atom.c*(_rfermi[0][n-1] - y[0]) + atom.a*(_rfermi[1][n-1] - y[1]));
    z *= atom.b;
    return z;
  }
}

double GetAtomicR(void) {
  return atom.rn;
}
