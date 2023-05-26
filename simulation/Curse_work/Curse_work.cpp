#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <array>
#include <algorithm>
#include <numeric>
#include <string>
#include <clocale>

using namespace std;

// Константы

const double     PI = 3.141592653589793;            // Число Пи
const double     TWO_PI = 2.0 * PI;                 // удвоенное число Пи
const double     E_CHARGE = 1.60217662e-19;         // заряд электрона [C]
const double     EV_TO_J = E_CHARGE;                // конвертационный фактор eV <-> Joule
const double     E_MASS = 9.10938356e-31;           // масса электрона [kg]
const double     AR_MASS = 6.63352090e-26;          // масса атома аргона [kg]
const double     MU_ARAR = AR_MASS / 2.0;           // приведенная масса двух атомов аргона [kg]
const double     K_BOLTZMANN = 1.38064852e-23;      // постоянная Больцмана [J/K]
const double     EPSILON0 = 8.85418781e-12;         // диэлектрическая проницаемость свободного пространства [F/m]

// параметры для симуляции

const int        N_G = 400;                         // количество точек в сетке
const int        N_T = 4000;                        // временные шаги в пределах периода радиочастоты RF
const double     FREQUENCY = 13.56e6;               // частота движения [Hz]
const double     VOLTAGE = 500.0;                   // амплитуда напряжения [V]
const double     L = 0.025;                         // межэлектродный зазор [m]
const double     PRESSURE = 100.0;                  // давление газа [Pa]
const double     TEMPERATURE = 300.0;               // температура фонового газа [K]
const double     WEIGHT = 7.0e4;                    // масса суперчастиц
const double     ELECTRODE_AREA = 1.0e-4;           // (фиктивная) площадь электрода [m^2]
const int        N_INIT = 1000;                     // количество начальных электронов и ионов

// дополнительные (производные) константы

const double     PERIOD = 1.0 / FREQUENCY;                              // длина периода радиочастоты (RF) [s]
const double     DT_E = PERIOD / static_cast<double>(N_T);              // временной шаг для електрона [s]
const int        N_SUB = 20;                                            // ионы движутся только в этих циклах (субциклирование)
const double     DT_I = N_SUB * DT_E;                                   // временной шаг для иона [s]
const double     DX = L / static_cast<double>(N_G - 1);                 // разделение пространственной сетки [m]
const double     INV_DX = 1.0 / DX;                                     // величина, обратная размеру пространственной сетки [1/m]
const double     GAS_DENSITY = PRESSURE / (K_BOLTZMANN * TEMPERATURE);  // фоновая плотность газа [1/m^3]
const double     OMEGA = TWO_PI * FREQUENCY;                            // угловая частота [rad/s]

// поперечное сечение электронов и ионов

const int        N_CS = 5;                          // общее количество процессов / поперечных сечений
const int        E_ELA = 0;                         // идентификатор процесса: электронный/эластичный
const int        E_EXC = 1;                         // идентификатор процесса: электрон/возбуждение
const int        E_ION = 2;                         // идентификатор процесса: электрон/ионизация
const int        I_ISO = 3;                         // идентификатор процесса: ионный/эластичный/изотропный
const int        I_BACK = 4;                        // идентификатор процесса: ионный/эластичный/обратное рассеяние
const double     E_EXC_TH = 11.5;                   // порог возбуждения электронным ударом [eV]
const double     E_ION_TH = 15.8;                   // порог ионизации при электронном ударе [eV]
const int        CS_RANGES = 1000000;               // количество вхождений в массивах поперечных сечений
const double     DE_CS = 0.001;                     // разделение энергии в массивах поперечного сечения [eV]
using cross_section = array<float, CS_RANGES>;      // массив поперечных сечений

cross_section    sigma[N_CS];                               // набор массивов поперечных сечений
cross_section    sigma_tot_e;                               // полное макроскопическое поперечное сечение электронов
cross_section    sigma_tot_i;                               // полное макроскопическое поперечное сечение ионов

// координаты частиц

size_t     N_e = 0;                                         // количество электронов
size_t     N_i = 0;                                         // количество ионов

vector<double>  x_e, vx_e{}, vy_e{}, vz_e{};                // координаты электронов (одна пространственная, три компоненты скорости)
vector<double>  x_i, vx_i, vy_i, vz_i;                      // координаты ионов (одна пространственная, три компоненты скорости)

using xvector = array<double, N_G>;                         // массив для величин, определенных в точках сетки
xvector  efield, pot;                                       // электрическое поле и потенциал
xvector  e_density, i_density;                              // плотности электронов и ионов
xvector  cumul_e_density, cumul_i_density;                  // совокупные плотности

using Ullong = unsigned long long int;                      // компактное имя для 64-разрядного целого числа без знака
Ullong   N_e_abs_pow = 0;                                   // счетчик электронов, поглощенных на питаемом электроде
Ullong   N_e_abs_gnd = 0;                                   // счетчик электронов, поглощенных на заземленном электроде
Ullong   N_i_abs_pow = 0;                                   // счетчик ионов, поглощенных на питаемом электроде
Ullong   N_i_abs_gnd = 0;                                   // счетчик ионов, поглощенных на заземленном электроде

// функция вероятности энергии электрона

const int    N_EEPF = 2000;                                 // количество энергетических ячеек в функции вероятности энергии электрона (EEPF)
const double DE_EEPF = 0.05;                                // резолюция EEPF [eV]
using eepf_vector = array<double, N_EEPF>;                  // массив для EEPF
eepf_vector  eepf{ 0.0 };                                   // интегрированная по времени EEPF в центре плазмы

// поток ионов - распределение энергии

const int    N_IFED = 200;                                  // количество энергетических ячеек в распределениях энергии потока ионов (Ion Flux-Energy Distributions - IFEDs)
const double DE_IFED = 1.0;                                 // разрешение IFEDs [eV]
using ifed_vector = array<int, N_IFED>;                     // массив для IFEDs
ifed_vector  ifed_pow{ 0 };                                 // IFED на питаемом электроде
ifed_vector  ifed_gnd{ 0 };                                 // IFED на заземленном электроде
double       mean_i_energy_pow;                             // средняя энергия ионов на питаемом электроде
double       mean_i_energy_gnd;                             // средняя энергия ионов на заземленном электроде

// пространственно-временные (XT) распределения

const int N_BIN = 20;                                       // количество временных шагов, привязанных к распределениям XT
const int N_XT = N_T / N_BIN;                               // количество пространственных ячеек для распределений XT
using xt_distr = array<double, N_G* N_XT>;                  // массив для распределений XT (десятичные числа)


xt_distr pot_xt = { 0.0 };                                  // XT распределение потенциала
xt_distr efield_xt = { 0.0 };                               // XT распределение электрического поля
xt_distr ne_xt = { 0.0 };                                   // XT распределение плотности электронов
xt_distr ni_xt = { 0.0 };                                   // XT распределение плотности ионов
xt_distr ue_xt = { 0.0 };                                   // XT распределение средней скорости электронов
xt_distr ui_xt = { 0.0 };                                   // XT распределение средней скорости ионов
xt_distr je_xt = { 0.0 };                                   // XT распределение плотности электронного тока
xt_distr ji_xt = { 0.0 };                                   // XT распределение плотности ионного тока
xt_distr powere_xt = { 0.0 };                               // XT распределение скорости возбуждения электронов (поглощения энергии)
xt_distr poweri_xt = { 0.0 };                               // XT распределение скорости ионного питания (поглощения энергии)
xt_distr meanee_xt = { 0.0 };                               // XT распределение средней энергии электронов
xt_distr meanei_xt = { 0.0 };                               // XT распределение средней энергии ионов
xt_distr counter_e_xt = { 0.0 };                            // XT счетчик свойств электронов
xt_distr counter_i_xt = { 0.0 };                            // XT счетчик свойств ионов
xt_distr ioniz_rate_xt = { 0.0 };                           // XT распределение скорости ионизации

double    mean_energy_accu_center = 0;                      // накопитель средней энергии электронов в центре зазора
Ullong    mean_energy_counter_center = 0;                   // счетчик средней энергии электронов в центре зазора
Ullong    N_e_coll = 0;                                     // счетчик электронных столкновений
Ullong    N_i_coll = 0;                                     // счетчик для ионных столкновений

double    Time;                                             // общее время моделирования (с начала моделирования)
int       cycle, no_of_cycles, cycles_done;                 // текущий цикл и общее количество циклов в прогоне, завершенные циклы (с начала моделирования)

int       arg1;                                             // для чтения аргументов командной строки
bool      measurement_mode;                                 // флаг, управляющий измерениями и сохранением данных

ofstream  datafile("conv.dat", ios_base::app);              // поток во внешний файл для сохранения данных конвергенции

//---------------------------------------------------------------------------------//
// Генератор C++ Mersenne Twister 19937                                            //
// R01 (MT gen) будет генерировать равномерное распределение по интервалу [0,1)    //
// RMB (MT gen) будет генерировать распределение Максвелла-Больцмана (атомов газа) //
//---------------------------------------------------------------------------------//

std::random_device rd{};
std::mt19937 MTgen(rd());
std::uniform_real_distribution<> R01(0.0, 1.0);
std::normal_distribution<> RMB(0.0, sqrt(K_BOLTZMANN* TEMPERATURE / AR_MASS));

//--------------------------------------------------------------------------------------//
//  поперечные сечения электронов: А. В. Фелпс и З. Л.Дж. Петрович, PSST 8 R21 (1999)   //
//--------------------------------------------------------------------------------------//

void set_electron_cross_sections_ar(void) {
    cout << ">> Состояние: настройка поперечных сечений e- / Ar" << endl;

    auto qmel = [](auto en) { 
        return 1e-20 * (fabs(6.0 / pow(1.0 + (en / 0.1) + pow(en / 0.6, 2.0), 3.3)
        - 1.1 * pow(en, 1.4) / (1.0 + pow(en / 15.0, 1.2)) / sqrt(1.0 + pow(en / 5.5, 2.5) + pow(en / 60.0, 4.1)))
        + 0.05 / pow(1.0 + en / 10.0, 2.0) + 0.01 * pow(en, 3.0) / (1.0 + pow(en / 12.0, 6.0))); 
    };

    auto qexc = [](const auto& en) { 
        if (en > E_EXC_TH) {
            return 1e-20 * (0.034 * pow(en - 11.5, 1.1) * (1.0 + pow(en / 15.0, 2.8)) / (1.0 + pow(en / 23.0, 5.5))
            + 0.023 * (en - 11.5) / pow(1.0 + en / 80.0, 1.9));
        }
        else { 
            return 0.0; 
        } 
    };

    auto qion = [](const auto& en) {
        if (en > E_ION_TH) {
            return 1e-20 * (970.0 * (en - 15.8) / pow(70.0 + en, 2.0) +
            0.06 * pow(en - 15.8, 2.0) * exp(-en / 9));
        }
        else { 
            return 0.0; 
        } 
    };

    vector<float> e(CS_RANGES);
    e[0] = DE_CS;
    generate(e.begin() + 1, e.end(), [i = 1]()mutable {
        return DE_CS * (i++);                                       // энергия электронов
    });                                            

    transform(e.begin(), e.end(), sigma[E_ELA].begin(), qmel);      // поперечное сечение для упругого столкновения e-/Ar
    transform(e.begin(), e.end(), sigma[E_EXC].begin(), qexc);      // поперечное сечение для возбуждения e-/Ar
    transform(e.begin(), e.end(), sigma[E_ION].begin(), qion);      // поперечное сечение для ионизации e-/Ar
}

//---------------------------------------------------------------------------//
//  поперечные сечения ионов: А. В. Фелпс, Дж. Аппл. Физика. 76, 747 (1994)  //
//---------------------------------------------------------------------------//

void set_ion_cross_sections_ar(void) {
    cout << ">> Состояние: настройка поперечных сечений Ar+ / Ar" << endl;
    auto qiso = [](const auto& e_lab) {
        return 2e-19 * pow(e_lab, -0.5) / (1.0 + e_lab) +
        3e-19 * e_lab / pow(1.0 + e_lab / 3.0, 2.0); 
    };

    auto qmom = [](const auto& e_lab) {
        return 1.15e-18 * pow(e_lab, -0.1) * pow(1.0 + 0.015 / e_lab, 0.6); 
    };

    auto qback = [&](const auto& x) { 
        return (qmom(x) - qiso(x)) / 2.0; 
    };

    vector<float> e(CS_RANGES);
    e[0] = 2.0 * DE_CS;
    generate(e.begin() + 1, e.end(), [i = 1]()mutable {return 2.0 * DE_CS * (i++); });   // энергия ионов в лабораторной системе отсчета

    transform(e.begin(), e.end(), sigma[I_ISO].begin(), qiso);     // поперечное сечение для Ar+/ Ar изотропной части упругого рассеяния
    transform(e.begin(), e.end(), sigma[I_BACK].begin(), qback);   // поперечное сечение для обратного упругого рассеяния Ar+ / Ar
}

//----------------------------------------------------------------------//
//  расчет полных поперечных сечений для электронов и ионов             //
//----------------------------------------------------------------------//

void calc_total_cross_sections(void) {

    for (size_t i{ 0 }; i < CS_RANGES; ++i) {
        sigma_tot_e[i] = (sigma[E_ELA][i] + sigma[E_EXC][i] + sigma[E_ION][i]) * GAS_DENSITY;   // полное макроскопическое поперечное сечение электронов
        sigma_tot_i[i] = (sigma[I_ISO][i] + sigma[I_BACK][i]) * GAS_DENSITY;                    // общее макроскопическое поперечное сечение ионов
    }
}

//----------------------------------------------------------------------//
//  испытание поперечных сечений для электронов и ионов                 //
//----------------------------------------------------------------------//

void test_cross_sections(void) {
    ofstream f("cross_sections.dat");                             // поперечные сечения, сохраненные в файле данных: cross_sections.data
    ostream_iterator<float> tofile(f, "\n");

    for (size_t i{ 0 }; i < CS_RANGES;++i) { f << i * DE_CS << endl; }
    for (const auto& v : sigma) {
        copy(v.begin(), v.end(), tofile);
    }
    f.close();
}

//---------------------------------------------------------------------//
// поиск верхнего предела частоты столкновений                         //
//---------------------------------------------------------------------//

double max_electron_coll_freq(void) {
    double e, v, nu, nu_max;
    nu_max = 0;
    for (size_t i{ 0 }; i < CS_RANGES; ++i) {
        e = i * DE_CS;
        v = sqrt(2.0 * e * EV_TO_J / E_MASS);
        nu = v * sigma_tot_e[i];
        if (nu > nu_max) { nu_max = nu; }
    }
    return nu_max;
}

double max_ion_coll_freq(void) {
    double e, g, nu, nu_max;
    nu_max = 0;
    for (size_t i{ 0 }; i < CS_RANGES; ++i) {
        e = i * DE_CS;
        g = sqrt(2.0 * e * EV_TO_J / MU_ARAR);
        nu = g * sigma_tot_i[i];
        if (nu > nu_max) nu_max = nu;
    }
    return nu_max;
}

//----------------------------------------------------------------------//
// инициализация моделирования путем размещения заданного количества    //
// электронов и ионов в случайных положениях между электродами          //
//----------------------------------------------------------------------//

void init(int nseed) {
    N_e = nseed;
    N_i = nseed;
    x_e.resize(nseed);
    x_i.resize(nseed);
    vx_e.resize(nseed, 0.0);                                      // компоненты начальной скорости электрона
    vy_e.resize(nseed, 0.0);
    vz_e.resize(nseed, 0.0);
    vx_i.resize(nseed, 0.0);                                      // компоненты начальной скорости иона
    vy_i.resize(nseed, 0.0);
    vz_i.resize(nseed, 0.0);
    generate(x_e.begin(), x_e.end(), [=]() {return L * R01(MTgen);});  // иницилизация начального случайного положения электрона
    generate(x_i.begin(), x_i.end(), [=]() {return L * R01(MTgen);});  // иницилизация начального случайного положения иона
}

//----------------------------------------------------------------------//
// e / Ar столкновение (приближение холодного газа)                     //
//----------------------------------------------------------------------//

void collision_electron(double xe, double& vxe, double& vye, double& vze, const int& eindex) {
    const double F1 = E_MASS / (E_MASS + AR_MASS);
    const double F2 = AR_MASS / (E_MASS + AR_MASS);
    double t0, t1, t2, rnd;
    double g, g2, gx, gy, gz, gx_old, gy_old, gz_old, wx, wy, wz, theta, phi;
    double chi, eta, chi2, eta2, sc, cc, se, ce, st, ct, sp, cp, energy, e_sc, e_ej;

    // вычисление относительной скорости перед столкновением и скорости центра масс
    gx = vxe;
    gy = vye;
    gz = vze;
    g = sqrt(gx * gx + gy * gy + gz * gz);
    wx = F1 * vxe;
    wy = F1 * vye;
    wz = F1 * vze;

    // нахождение угла Эйлера:

    if (gx == 0) { theta = 0.5 * PI; }
    else { theta = atan2(sqrt(gy * gy + gz * gz), gx); }
    if (gy == 0) {
        if (gz > 0) { phi = 0.5 * PI; }
        else { phi = -0.5 * PI; }
    }
    else { phi = atan2(gz, gy); }
    st = sin(theta);
    ct = cos(theta);
    sp = sin(phi);
    cp = cos(phi);

    // выбор типа столкновения на основе поперечных сечений
    // учитывание потерь энергии при неупругих столкновениях
    // генерация углов рассеяния и азимута
    // в случае ионизации оперирование с "новым" электроном

    t0 = sigma[E_ELA][eindex];
    t1 = t0 + sigma[E_EXC][eindex];
    t2 = t1 + sigma[E_ION][eindex];
    rnd = R01(MTgen);
    if (rnd < (t0 / t2)) {                          // упругое рассеяние
        chi = acos(1.0 - 2.0 * R01(MTgen));         // изотропное рассеяние
        eta = TWO_PI * R01(MTgen);                  // азимутальный угол
    }
    else if (rnd < (t1 / t2)) {                     // возбуждение
        energy = 0.5 * E_MASS * g * g;              // энергия электрона
        energy = fabs(energy - E_EXC_TH * EV_TO_J); // потеря энергии на возбуждение
        g = sqrt(2.0 * energy / E_MASS);            // относительная скорость после потери энергии
        chi = acos(1.0 - 2.0 * R01(MTgen));         // изотропное рассеяние
        eta = TWO_PI * R01(MTgen);                  // азимутальный угол
    }
    else {                                          // иницилизация
        energy = 0.5 * E_MASS * g * g;              // энергия электрона
        energy = fabs(energy - E_ION_TH * EV_TO_J); // потеря энергии на ионизацию
        e_ej = 10.0 * tan(R01(MTgen) * atan(energy / EV_TO_J / 20.0)) * EV_TO_J;    // энергия испускаемого электрона
        e_sc = fabs(energy - e_ej);                 // энергия входящего электрона после столкновения
        g = sqrt(2.0 * e_sc / E_MASS);              // относительная скорость входящего (исходного) электрона
        g2 = sqrt(2.0 * e_ej / E_MASS);             // относительная скорость испущенного (нового) электрона
        chi = acos(sqrt(e_sc / energy));            // угол рассеяния для входящего электрона
        chi2 = acos(sqrt(e_ej / energy));           // угол рассеяния для испускаемых электронов
        eta = TWO_PI * R01(MTgen);                  // азимутальный угол для входящего электрона
        eta2 = eta + PI;                            // азимутальный угол для испускаемого электрона
        sc = sin(chi2);
        cc = cos(chi2);
        se = sin(eta2);
        ce = cos(eta2);
        gx = g2 * (ct * cc - st * sc * ce);
        gy = g2 * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
        gz = g2 * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
        N_e++;                                        // добавление нового електрона
        x_e.push_back(xe);
        N_i++;                                        // добавление нового иона
        x_i.push_back(xe);
        vx_i.push_back(RMB(MTgen));                   // скорость выбирается из фонового теплового распределения
        vy_i.push_back(RMB(MTgen));
        vz_i.push_back(RMB(MTgen));
    }
    // сохранение старых данных

    gx_old = wx + F2 * gx;
    gy_old = wy + F2 * gy;
    gz_old = wz + F2 * gz;
    
    // рассеяние первичного электрона

    sc = sin(chi);
    cc = cos(chi);
    se = sin(eta);
    ce = cos(eta);

    // вычисление новой относительной скорости:

    gx = g * (ct * cc - st * sc * ce);
    gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);

    // скорость электрона после столкновения

    vxe = wx + F2 * gx;
    vye = wy + F2 * gy;
    vze = wz + F2 * gz;

    if (rnd < (t0 / t2)) {
    }
    else if (rnd < (t1 / t2)) {
    }
    else {
        vx_e.push_back(wx + F2 * gx);
        vy_e.push_back(wy + F2 * gy);
        vz_e.push_back(wz + F2 * gz);
    }
}
//----------------------------------------------------------------------//
// Ar+ / Ar столкновения                                                //
//----------------------------------------------------------------------//

void collision_ion(double& vx_1, double& vy_1, double& vz_1,
    double& vx_2, double& vy_2, double& vz_2, const int& e_index) {
    double   g, gx, gy, gz, wx, wy, wz, rnd;
    double   theta, phi, chi, eta, st, ct, sp, cp, sc, cc, se, ce, t1, t2;

    // вычисление относительной скорости до столкновения
    // случайный целевой атом Максвелла уже выбран 
    // (компоненты скорости vx_2,vy_2,vz_2 целевого атома поступают вместе с вызовом)


    gx = vx_1 - vx_2;
    gy = vy_1 - vy_2;
    gz = vz_1 - vz_2;
    g = sqrt(gx * gx + gy * gy + gz * gz);
    wx = 0.5 * (vx_1 + vx_2);
    wy = 0.5 * (vy_1 + vy_2);
    wz = 0.5 * (vz_1 + vz_2);

    // поиск угла Эйлера:

    if (gx == 0) { theta = 0.5 * PI; }
    else { theta = atan2(sqrt(gy * gy + gz * gz), gx); }
    if (gy == 0) {
        if (gz > 0) { phi = 0.5 * PI; }
        else { phi = -0.5 * PI; }
    }
    else { phi = atan2(gz, gy); }


    // определение типа столкновения на основе поперечных сечений и генерация угола рассеяния

    t1 = sigma[I_ISO][e_index];
    t2 = t1 + sigma[I_BACK][e_index];
    rnd = R01(MTgen);
    if (rnd < (t1 / t2)) {                          // изотропное рассеяние
        chi = acos(1.0 - 2.0 * R01(MTgen));         // угол изотропного рассеяния
    }
    else {                                          // обратное рассеяние
        chi = PI;                                   // угол обратного рассеяния
    }
    eta = TWO_PI * R01(MTgen);                      // азимутальный угол
    sc = sin(chi);
    cc = cos(chi);
    se = sin(eta);
    ce = cos(eta);
    st = sin(theta);
    ct = cos(theta);
    sp = sin(phi);
    cp = cos(phi);

    // вычисление новой относительной скорости:

    gx = g * (ct * cc - st * sc * ce);
    gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);

    // скорость иона после столкновения

    vx_1 = wx + 0.5 * gx;
    vy_1 = wy + 0.5 * gy;
    vz_1 = wz + 0.5 * gz;
}

//----------------------------------------------------------------------//
// решение уравнение Пуассона (алгоритм Томаса)                         //
//----------------------------------------------------------------------//

void solve_Poisson(const xvector& rho1, const double& tt) {
    const double A = 1.0;
    const double B = -2.0;
    const double C = 1.0;
    const double S = 1.0 / (2.0 * DX);
    const double ALPHA = -DX * DX / EPSILON0;
    xvector g{}, w{}, f{};
    size_t  i;

    // приложение потенциала к электродам - граничные условия

    pot.front() = VOLTAGE * cos(OMEGA * tt);    // потенциал на питаемом электроде
    pot.back() = 0.0;                           // потенциал на заземленном электроде

    // решение уравнения Пуассона

    for (i = 1; i <= N_G - 2; ++i) f[i] = ALPHA * rho1[i];
    f[1] -= pot.front();
    f[N_G - 2] -= pot.back();
    w[1] = C / B;
    g[1] = f[1] / B;
    for (i = 2; i <= N_G - 2; ++i) {
        w[i] = C / (B - A * w[i - 1]);
        g[i] = (f[i] - A * g[i - 1]) / (B - A * w[i - 1]);
    }
    pot[N_G - 2] = g[N_G - 2];
    for (i = N_G - 3; i > 0; --i) pot[i] = g[i] - w[i] * pot[i + 1];    // потенциал в точках сетки между электродами

    // вычисление электрического поля

    for (i = 1; i <= N_G - 2; ++i) efield[i] = (pot[i - 1] - pot[i + 1]) * S;                       // электрическое поле в точках сетки между электродами
    efield.front() = (pot[0] - pot[1]) * INV_DX - rho1.front() * DX / (2.0 * EPSILON0);             // питаемый электрод
    efield.back() = (pot[N_G - 2] - pot[N_G - 1]) * INV_DX + rho1.back() * DX / (2.0 * EPSILON0);   // заземленный электрод
}

//---------------------------------------------------------------------//
// моделирование одного радиочастотного цикла                          //
//---------------------------------------------------------------------//

void do_one_cycle(void) {
    const double DV = ELECTRODE_AREA * DX;
    const double FACTOR_W = WEIGHT / DV;
    const double FACTOR_E = DT_E / E_MASS * E_CHARGE;
    const double FACTOR_I = DT_I / AR_MASS * E_CHARGE;
    const double MIN_X = 0.45 * L;                      // минимальное положение для сбора EEPF
    const double MAX_X = 0.55 * L;                      // максимальное положение для сбора EEPF
    size_t       k, t, p, energy_index;
    double       rmod, rint, g, g_sqr, gx, gy, gz, vx_a, vy_a, vz_a, e_x, energy, nu, p_coll, v_sqr, velocity;
    double       mean_v, rate;
    bool         out;
    xvector      rho{};
    size_t       t_index;

    for (t = 0; t < N_T; t++) {                         // период RF делится на N_T равных временных интервалов (временной шаг DT_E)
        Time += DT_E;                                   // обносвление общего времени моделирования
        t_index = t / N_BIN;                            // индекс для распределений XT

        // Шаг 1: вычисление плотности в точках сетки

        fill(e_density.begin(), e_density.end(), 0.0);          // электронная плотность - вычисляется на каждом временном шаге
        for (k = 0; k < N_e; ++k) {
            rmod = modf(x_e[k] * INV_DX, &rint);
            p = static_cast<int>(rint);
            e_density[p] += (1.0 - rmod) * FACTOR_W;
            e_density[p + 1] += rmod * FACTOR_W;
        }
        e_density.front() *= 2.0;
        e_density.back() *= 2.0;
        transform(cumul_e_density.begin(), cumul_e_density.end(), e_density.begin(), cumul_e_density.begin(), [](auto x, auto y) {return x + y;});

        if ((t % N_SUB) == 0) {                                 // плотность ионов - вычисляется на каждом N_SUB-ом временном шаге (субциклирование)
            fill(i_density.begin(), i_density.end(), 0.0);
            for (k = 0; k < N_i; ++k) {
                rmod = modf(x_i[k] * INV_DX, &rint);
                p = static_cast<int>(rint);
                i_density[p] += (1.0 - rmod) * FACTOR_W;
                i_density[p + 1] += rmod * FACTOR_W;
            }
            i_density.front() *= 2.0;
            i_density.back() *= 2.0;
        }
        transform(cumul_i_density.begin(), cumul_i_density.end(), i_density.begin(), cumul_i_density.begin(), [](auto x, auto y) {return x + y;});

        // Шаг 2: решение уравнения Пуассона

        // получение плотности заряда
        transform(i_density.begin(), i_density.end(), e_density.begin(), rho.begin(), [](auto x, auto y) {return E_CHARGE * (x - y);});
        solve_Poisson(rho, Time);                               // вычисление потенциала и электрического поля

        // Шаги 3, 4: перемещение частицы в соответствии с электрическим полем, интерполированным в положение частиц

        for (k = 0; k < N_e; k++) {                             // перемещение всех электронов на каждом временном шаге
            rmod = modf(x_e[k] * INV_DX, &rint);
            p = static_cast<int>(rint);
            e_x = (1.0 - rmod) * efield[p] + rmod * efield[p + 1];

            if (measurement_mode) {

                // измерение: "x" и "v" необходимы одновременно, т.е. старое "x" и среднее "v"

                mean_v = vx_e[k] - 0.5 * e_x * FACTOR_E;
                counter_e_xt[p * N_XT + t_index] += (1.0 - rmod);
                counter_e_xt[(p + 1) * N_XT + t_index] += rmod;
                ue_xt[p * N_XT + t_index] += (1.0 - rmod) * mean_v;
                ue_xt[(p + 1) * N_XT + t_index] += rmod * mean_v;
                v_sqr = mean_v * mean_v + vy_e[k] * vy_e[k] + vz_e[k] * vz_e[k];
                energy = 0.5 * E_MASS * v_sqr / EV_TO_J;
                meanee_xt[p * N_XT + t_index] += (1.0 - rmod) * energy;
                meanee_xt[(p + 1) * N_XT + t_index] += rmod * energy;
                energy_index = min(static_cast<int>(energy / DE_CS + 0.5), CS_RANGES - 1);
                velocity = sqrt(v_sqr);
                rate = sigma[E_ION][energy_index] * velocity * DT_E * GAS_DENSITY;
                ioniz_rate_xt[p * N_XT + t_index] += (1.0 - rmod) * rate;
                ioniz_rate_xt[(p + 1) * N_XT + t_index] += rmod * rate;

                // измерение EEPF в центре

                if ((MIN_X < x_e[k]) && (x_e[k] < MAX_X)) {
                    energy_index = static_cast<int>(energy / DE_EEPF);
                    if (energy_index < N_EEPF) { eepf[energy_index] += 1.0; }
                    mean_energy_accu_center += energy;
                    mean_energy_counter_center++;
                }
            }

            // обновление скорости и положения

            vx_e[k] -= e_x * FACTOR_E;
            x_e[k] += vx_e[k] * DT_E;
        }

        if ((t % N_SUB) == 0) {                    // перемещение всех ионов за каждый N_SUB-й временной шаг (субциклирование)
            for (k = 0; k < N_i; k++) {
                rmod = modf(x_i[k] * INV_DX, &rint);
                p = static_cast<int>(rint);
                e_x = (1.0 - rmod) * efield[p] + rmod * efield[p + 1];

                if (measurement_mode) {

                    // измерение: "x" и "v" необходимы одновременно, т.е. старое "x" и среднее "v"

                    mean_v = vx_i[k] + 0.5 * e_x * FACTOR_I;
                    counter_i_xt[p * N_XT + t_index] += (1.0 - rmod);
                    counter_i_xt[(p + 1) * N_XT + t_index] += rmod;
                    ui_xt[p * N_XT + t_index] += (1.0 - rmod) * mean_v;
                    ui_xt[(p + 1) * N_XT + t_index] += rmod * mean_v;
                    v_sqr = mean_v * mean_v + vy_i[k] * vy_i[k] + vz_i[k] * vz_i[k];
                    energy = 0.5 * AR_MASS * v_sqr / EV_TO_J;
                    meanei_xt[p * N_XT + t_index] += (1.0 - rmod) * energy;
                    meanei_xt[(p + 1) * N_XT + t_index] += rmod * energy;
                }

                // обновление скорости и положения

                vx_i[k] += e_x * FACTOR_I;
                x_i[k] += vx_i[k] * DT_I;
            }
        }


        // Шаг 5: проверка границ
        k = 0;
        while (k < x_e.size()) {                            // проверка границ для всех электронов на каждом временном шаге
            out = false;
            if (x_e[k] < 0) { N_e_abs_pow++; out = true; }  // электрон отсутствует на питаемом электроде
            if (x_e[k] > L) { N_e_abs_gnd++; out = true; }  // электрон находится вне заземленного электрода
            if (out) {                                      // удаление электрона, если он отсутствует
                x_e[k] = x_e.back(); x_e.pop_back();
                vx_e[k] = vx_e.back(); vx_e.pop_back();
                vy_e[k] = vy_e.back(); vy_e.pop_back();
                vz_e[k] = vz_e.back(); vz_e.pop_back();
                N_e--;
            }
            else k++;
        }

        if ((t % N_SUB) == 0) {                             // проверка границ для всех ионов на каждом N_SUB-ом временном шаге (субциклирование)
            k = 0;
            while (k < x_i.size()) {
                out = false;
                if (x_i[k] < 0) {                           // выход иона на питаемый электрод
                    N_i_abs_pow++;
                    out = true;
                    v_sqr = vx_i[k] * vx_i[k] + vy_i[k] * vy_i[k] + vz_i[k] * vz_i[k];
                    energy = 0.5 * AR_MASS * v_sqr / EV_TO_J;
                    energy_index = static_cast<int>(energy / DE_IFED);
                    if (energy_index < N_IFED) { ifed_pow[energy_index]++; }    // сохраните IFED на питаемом электроде
                }
                if (x_i[k] > L) {                           // выход иона на заземленный электрод
                    N_i_abs_gnd++;
                    out = true;
                    v_sqr = vx_i[k] * vx_i[k] + vy_i[k] * vy_i[k] + vz_i[k] * vz_i[k];
                    energy = 0.5 * AR_MASS * v_sqr / EV_TO_J;
                    energy_index = static_cast<int>(energy / DE_IFED);
                    if (energy_index < N_IFED) { ifed_gnd[energy_index]++; }    // сохранение IFED на заземленном электроде
                }
                if (out) {                                  // удалите ион, если он отсутствует
                    x_i[k] = x_i.back(); x_i.pop_back();
                    vx_i[k] = vx_i.back(); vx_i.pop_back();
                    vy_i[k] = vy_i.back(); vy_i.pop_back();
                    vz_i[k] = vz_i.back(); vz_i.pop_back();
                    N_i--;
                }
                else k++;
            }
        }

        // Шаг 6: столкновение

        for (k = 0; k < N_e; ++k) {                     // проверка возникновения столкновения для всех электронов на каждом временном шаге
            v_sqr = vx_e[k] * vx_e[k] + vy_e[k] * vy_e[k] + vz_e[k] * vz_e[k];
            velocity = sqrt(v_sqr);
            energy = 0.5 * E_MASS * v_sqr / EV_TO_J;
            energy_index = min(static_cast<int>(energy / DE_CS + 0.5), CS_RANGES - 1);
            nu = sigma_tot_e[energy_index] * velocity;
            p_coll = 1 - exp(-nu * DT_E);               // вероятность столкновения электронов
            if (R01(MTgen) < p_coll) {                  // столкновение электронов
                collision_electron(x_e[k], vx_e[k], vy_e[k], vz_e[k], energy_index);
                N_e_coll++;
            }
        }

        if ((t % N_SUB) == 0) {                         // проверка возникновения столкновения для всех ионов на каждом N_SUB-ом временном шаге (субциклирование)
            for (k = 0; k < N_i; ++k) {
                vx_a = RMB(MTgen);                      // выбор компонентов скорости случайных атомов газа
                vy_a = RMB(MTgen);
                vz_a = RMB(MTgen);
                gx = vx_i[k] - vx_a;                    // вычисление относительной скорости партнеров по столкновению
                gy = vy_i[k] - vy_a;
                gz = vz_i[k] - vz_a;
                g_sqr = gx * gx + gy * gy + gz * gz;
                g = sqrt(g_sqr);
                energy = 0.5 * MU_ARAR * g_sqr / EV_TO_J;
                energy_index = min(static_cast<int>(energy / DE_CS + 0.5), CS_RANGES - 1);
                nu = sigma_tot_i[energy_index] * g;
                p_coll = 1 - exp(-nu * DT_I);           // вероятность столкновения ионов
                if (R01(MTgen) < p_coll) {              // столкновение ионов
                    collision_ion(vx_i[k], vy_i[k], vz_i[k], vx_a, vy_a, vz_a, energy_index);
                    N_i_coll++;
                }
            }
        }

        if (measurement_mode) {

            // сбор данных из сетки:

            for (p = 0; p < N_G; p++) {
                pot_xt[p * N_XT + t_index] += pot[p];
                efield_xt[p * N_XT + t_index] += efield[p];
                ne_xt[p * N_XT + t_index] += e_density[p];
                ni_xt[p * N_XT + t_index] += i_density[p];
            }
        }

        if ((t % 1000) == 0) {
            cout << " c = " << setw(8) << cycle << "  t = " << setw(8) << t << "  #e = " << setw(8) << N_e << "  #i = " << setw(8) << N_i << endl;
        }
    }
    datafile << cycle << "\t" << N_e << "\t" << N_i << endl;
}

//---------------------------------------------------------------------//
// сохранение и загрузка координат частиц                              //
//---------------------------------------------------------------------//

void save_particle_data() {
    ofstream f("picdata.bin", ios::binary);

    f.write(reinterpret_cast<char*>(&Time), sizeof(double));
    f.write(reinterpret_cast<char*>(&cycles_done), sizeof(int));
    f.write(reinterpret_cast<char*>(&N_e), sizeof(int));
    f.write(reinterpret_cast<char*>(&x_e[0]), N_e * sizeof(double));
    f.write(reinterpret_cast<char*>(&vx_e[0]), N_e * sizeof(double));
    f.write(reinterpret_cast<char*>(&vy_e[0]), N_e * sizeof(double));
    f.write(reinterpret_cast<char*>(&vz_e[0]), N_e * sizeof(double));
    f.write(reinterpret_cast<char*>(&N_i), sizeof(int));
    f.write(reinterpret_cast<char*>(&x_i[0]), N_i * sizeof(double));
    f.write(reinterpret_cast<char*>(&vx_i[0]), N_i * sizeof(double));
    f.write(reinterpret_cast<char*>(&vy_i[0]), N_i * sizeof(double));
    f.write(reinterpret_cast<char*>(&vz_i[0]), N_i * sizeof(double));
    f.close();

    cout << ">> Состояние: сохранение данных : " << N_e << " электронов " << N_i << " ионов, ";
    if (cycles_done % 10 == 1)
        cout << cycles_done << " цикл выполнен, времени прошло " << scientific << Time << " [s]" << endl;
    else if (cycles_done % 10 >= 2 && cycles_done % 10 <= 4)
        cout << cycles_done << " цикла выполнена, времени прошло " << scientific << Time << " [s]" << endl;
    else
        cout << cycles_done << " циклов выполнено, времени прошло " << scientific << Time << " [s]" << endl;
}

//---------------------------------------------------------------------//
// загрузка координат частиц                                           //
//---------------------------------------------------------------------//

void load_particle_data() {
    ifstream f("picdata.bin", std::ios::binary);
    if (f.fail()) { cout << ">> ОШИБКА: файл данных о частицах не найден, попробуйте запустить начальный цикл, используя аргумент '0'" << endl; exit(0); }
    f.read(reinterpret_cast<char*>(&Time), sizeof(double));
    f.read(reinterpret_cast<char*>(&cycles_done), sizeof(int));
    f.read(reinterpret_cast<char*>(&N_e), sizeof(int));
    x_e.resize(N_e);
    vx_e.resize(N_e);
    vy_e.resize(N_e);
    vz_e.resize(N_e);
    f.read(reinterpret_cast<char*>(&x_e[0]), N_e * sizeof(double));
    f.read(reinterpret_cast<char*>(&vx_e[0]), N_e * sizeof(double));
    f.read(reinterpret_cast<char*>(&vy_e[0]), N_e * sizeof(double));
    f.read(reinterpret_cast<char*>(&vz_e[0]), N_e * sizeof(double));
    f.read(reinterpret_cast<char*>(&N_i), sizeof(int));
    x_i.resize(N_i);
    vx_i.resize(N_i);
    vy_i.resize(N_i);
    vz_i.resize(N_i);
    f.read(reinterpret_cast<char*>(&x_i[0]), N_i * sizeof(double));
    f.read(reinterpret_cast<char*>(&vx_i[0]), N_i * sizeof(double));
    f.read(reinterpret_cast<char*>(&vy_i[0]), N_i * sizeof(double));
    f.read(reinterpret_cast<char*>(&vz_i[0]), N_i * sizeof(double));
    f.close();

    cout << ">> Состояние: загрузка данных: " << N_e << " электронов " << N_i << " ионов, ";
    if (cycles_done%10 == 1)
        cout << cycles_done << " завершенный цикл, времени прошло " << scientific << Time << " [s]" << endl;
    else if (cycles_done % 10 >= 2 && cycles_done % 10 <= 4)
        cout << cycles_done << " завершенных цикла, времени прошло " << scientific << Time << " [s]" << endl;
    else
        cout << cycles_done << " завершенных циклов, времени прошло " << scientific << Time << " [s]" << endl;
}

//---------------------------------------------------------------------//
// сохранение данных о плотности                                       //
//---------------------------------------------------------------------//

void save_density(void) {
    ofstream f("density.dat");
    f << setprecision(12) << fixed << scientific;

    auto c = 1.0 / static_cast<double>(no_of_cycles) / static_cast<double>(N_T);
    for (size_t i{ 0 }; i < N_G;++i) {
        f << i * DX << "\t" << cumul_e_density[i] * c << "\t" << cumul_i_density[i] * c << endl;
    }
    f.close();
}

//---------------------------------------------------------------------//
// Сохранение EEPF данных                                              //
//---------------------------------------------------------------------//

void save_eepf(void) {
    ofstream f("eepf.dat");
    auto h = accumulate(eepf.begin(), eepf.end(), 0.0);
    h *= DE_EEPF;
    f << scientific;
    double energy{};
    for (size_t i{ 0 }; i < N_EEPF;++i) {
        energy = (i + 0.5) * DE_EEPF;
        f << energy << "\t" << eepf[i] / h / sqrt(energy) << endl;
    }
    f.close();
}

//---------------------------------------------------------------------//
// Сохранение IFED данных                                              //
//---------------------------------------------------------------------//

void save_ifed(void) {
    double p, g, energy;
    ofstream f("ifed.dat");
    f << scientific;
    double h_pow = accumulate(ifed_pow.begin(), ifed_pow.end(), 0.0);
    double h_gnd = accumulate(ifed_gnd.begin(), ifed_gnd.end(), 0.0);
    h_pow *= DE_IFED;
    h_gnd *= DE_IFED;
    mean_i_energy_pow = 0.0;
    mean_i_energy_gnd = 0.0;
    for (size_t i{ 0 }; i < N_IFED;++i) {
        energy = (i + 0.5) * DE_IFED;
        p = static_cast<double>(ifed_pow[i]) / h_pow;
        g = static_cast<double>(ifed_gnd[i]) / h_gnd;
        f << energy << "\t" << p << "\t" << g << endl;
        mean_i_energy_pow += energy * p;
        mean_i_energy_gnd += energy * g;
    }
    f.close();
}

//--------------------------------------------------------------------//
// Сохранение XT данных                                               //
//--------------------------------------------------------------------//

void save_xt_1(xt_distr& distr, string fname) {
    ofstream f(fname);
    ostream_iterator<double> tof(f, " ");
    auto it = distr.begin();

    f << setprecision(8) << fixed << scientific;
    for (size_t i{ 0 };i < N_G;++i) {
        copy_n(it, N_XT, tof);
        advance(it, N_XT);
        f << endl;
    }
    f.close();
}



void norm_all_xt(void) {
    // нормализование всех данных XT

    double f1 = static_cast<double>(N_XT) / static_cast<double>(no_of_cycles * N_T);
    double f2 = WEIGHT / (ELECTRODE_AREA * DX) / (no_of_cycles * (PERIOD / static_cast<double>(N_XT)));

    transform(pot_xt.begin(), pot_xt.end(), pot_xt.begin(), [=](auto y) {return f1 * y;});
    transform(efield_xt.begin(), efield_xt.end(), efield_xt.begin(), [=](auto y) {return f1 * y;});
    transform(ne_xt.begin(), ne_xt.end(), ne_xt.begin(), [=](auto y) {return f1 * y;});
    transform(ni_xt.begin(), ni_xt.end(), ni_xt.begin(), [=](auto y) {return f1 * y;});

    transform(ue_xt.begin(), ue_xt.end(), counter_e_xt.begin(), ue_xt.begin(), [](auto x, auto y) {
        if (y > 0) { 
            return x / y; } 
        else { 
            return 0.0; 
        }
        }
    );
    transform(ue_xt.begin(), ue_xt.end(), ne_xt.begin(), je_xt.begin(), [=](auto x, auto y) {
        return -x * y * E_CHARGE;
        }
    );
    transform(meanee_xt.begin(), meanee_xt.end(), counter_e_xt.begin(), meanee_xt.begin(), [](auto x, auto y) {
        if (y > 0) { 
            return x / y; } 
        else { 
            return 0.0; 
        }
        }
    );
    transform(ioniz_rate_xt.begin(), ioniz_rate_xt.end(), counter_e_xt.begin(), ioniz_rate_xt.begin(), [=](auto x, auto y) {
        if (y > 0) { 
            return x * f2; } 
        else { 
            return 0.0; 
        }
        }
    );

    transform(ui_xt.begin(), ui_xt.end(), counter_i_xt.begin(), ui_xt.begin(), [](auto x, auto y) {
        if (y > 0) { 
            return x / y; 
        } else { 
            return 0.0; 
        }
        }
    );
    transform(ui_xt.begin(), ui_xt.end(), ni_xt.begin(), ji_xt.begin(), [=](auto x, auto y) {
        return x * y * E_CHARGE;
        }
    );
    transform(meanei_xt.begin(), meanei_xt.end(), counter_i_xt.begin(), meanei_xt.begin(), [](auto x, auto y) {
        if (y > 0) { 
            return x / y; } 
        else { 
            return 0.0; 
        }
        }
    );

    transform(je_xt.begin(), je_xt.end(), efield_xt.begin(), powere_xt.begin(), [=](auto x, auto y) {
        return x * y;
        }
    );
    transform(ji_xt.begin(), ji_xt.end(), efield_xt.begin(), poweri_xt.begin(), [=](auto x, auto y) {
        return x * y;
        }
    );
}


void save_all_xt(void) {

    save_xt_1(pot_xt, "pot_xt.dat");
    save_xt_1(efield_xt, "efield_xt.dat");
    save_xt_1(ne_xt, "ne_xt.dat");
    save_xt_1(ni_xt, "ni_xt.dat");
    save_xt_1(je_xt, "je_xt.dat");
    save_xt_1(ji_xt, "ji_xt.dat");
    save_xt_1(powere_xt, "powere_xt.dat");
    save_xt_1(poweri_xt, "poweri_xt.dat");
    save_xt_1(meanee_xt, "meanee_xt.dat");
    save_xt_1(meanei_xt, "meanei_xt.dat");
    save_xt_1(ioniz_rate_xt, "ioniz_xt.dat");
}

//---------------------------------------------------------------------//
// отчет о моделировании, включающий условия стабильности и точности   //
//---------------------------------------------------------------------//

void check_and_save_info(void) {
    ofstream f("info.txt");
    string line(80, '-');
    f << setprecision(4) << fixed << scientific;

    double density = cumul_e_density[N_G / 2] / static_cast<double>(no_of_cycles) / static_cast<double>(N_T);   // e плотность @ center
    double plas_freq = E_CHARGE * sqrt(density / EPSILON0 / E_MASS);                                            // e плазменная частота @ center
    double meane = mean_energy_accu_center / static_cast<double>(mean_energy_counter_center);                   // e средняя энергия @ center
    double kT = 2.0 * meane * EV_TO_J / 3.0;                                                                    // k*T_e @ center (approximate)
    double debye_length = sqrt(EPSILON0 * kT / density) / E_CHARGE;                                             // e Дебаевский радиус @ center
    double sim_time = static_cast<double>(no_of_cycles) / FREQUENCY;                                            // время моделирования
    double ecoll_freq = static_cast<double>(N_e_coll) / sim_time / static_cast<double>(N_e);                    // частота столкновений электронов
    double icoll_freq = static_cast<double>(N_i_coll) / sim_time / static_cast<double>(N_i);                    // частота столкновений ионов

    f << "############################## Отчёт симуляции ################################" << endl;
    f << "Параметры моделирования:" << endl;
    f << "Расстояние между зазорами                     = " << L << " [m]" << endl;
    f << "Кол-во разделений сетки                       = " << N_G << endl;
    f << "Частота                                       = " << FREQUENCY << " [Hz]" << endl;
    f << "Кол-во временных шагов / период               = " << N_T << endl;
    f << "Кол-во временных шагов электрона/иона         = " << N_SUB << endl;
    f << "Амплитуда напряжения                          = " << VOLTAGE << " [V]" << endl;
    f << "Давление (Ar)                                 = " << PRESSURE << " [Pa]" << endl;
    f << "Температура                                   = " << TEMPERATURE << " [K]" << endl;
    f << "Масса супер-частицы                           = " << WEIGHT << endl;
    f << "Кол-во циклов моделирования в этом прогоне    = " << no_of_cycles << endl;
    f << line << endl;
    f << "Характеристики плазмы:" << endl;
    f << "Электронная плотность @ center                = " << density << " [m^{-3}]" << endl;
    f << "Плазменная частота @ center                   = " << plas_freq << " [rad/s]" << endl;
    f << "Дебаевский радиус @ center                    = " << debye_length << " [m]" << endl;
    f << "Частота столкновения электронов               = " << ecoll_freq << " [1/s]" << endl;
    f << "Частота столкновения ионов                    = " << icoll_freq << " [1/s]" << endl;
    f << line << endl;
    f << "Условия стабильности и точности:" << endl;
    auto conditions_OK = true;
    auto c = plas_freq * DT_E;
    f << "Плазменная частота @ center * DT_e            = " << c << " (OK, если менее 0.20)" << endl;
    if (c > 0.2) { conditions_OK = false; }
    c = DX / debye_length;
    f << "DX / Дебаевский радиус @ center               = " << c << " (OK, если менее 1.00)" << endl;
    if (c > 1.0) { conditions_OK = false; }
    c = max_electron_coll_freq() * DT_E;
    f << "Макс. частота столкновений электронов * DT_E  = " << c << " (OK, если менее 0.05)" << endl;
    if (c > 0.05) { conditions_OK = false; }
    c = max_ion_coll_freq() * DT_I;
    f << "Макс. частота столкновений ионов * DT_I       = " << c << " (OK, если менее 0.05)" << endl;
    if (c > 0.05) { conditions_OK = false; }
    if (conditions_OK == false) {
        f << line << endl;
        f << "** НАРУШЕНЫ УСЛОВИЯ СТАБИЛЬНОСТИ И ТОЧНОСТИ - УТОЧНИТЕ НАСТРОЙКИ МОДЕЛИРОВАНИЯ! **" << endl;
        f << line << endl;
        f.close();
        f << ">> Состояние: ОШИБКА: НАРУШЕНЫ УСЛОВИЯ СТАБИЛЬНОСТИ И ТОЧНОСТИ! " << endl;
        f << ">> Состояние: Для получения подробной информации смотрите 'info.txt' и уточняйте настройки моделирования!" << endl;
    }
    else {
        // вычисляет максимальную энергию, для которой выполняется условие Куранта:
        double v_max = DX / DT_E;
        double e_max = 0.5 * E_MASS * v_max * v_max / EV_TO_J;
        f << "Макс. e- энергия для условий CFL              = " << e_max << endl;
        f << "Проверьте EEPF, чтобы убедиться, что CFL выполняется для большинства электронов!" << endl;
        f << line << endl;

        // сохранение следующих данных выполняется здесь, так как некоторые из дальнейших строк нуждаются в новых данных,
        // которые вычисляются / нормализуются в этих функциях

        cout << ">> Состояние: сохранение диагностических данных" << endl;
        save_density();
        save_eepf();
        save_ifed();
        norm_all_xt();
        save_all_xt();
        f << "Характеристики частиц на электродах:" << endl;
        f << "Поток ионов на питаемом электроде                 = " << N_i_abs_pow * WEIGHT / ELECTRODE_AREA / (no_of_cycles * PERIOD) << " [m^{-2} s^{-1}]" << endl;
        f << "Поток ионов на заземленном электроде              = " << N_i_abs_gnd * WEIGHT / ELECTRODE_AREA / (no_of_cycles * PERIOD) << " [m^{-2} s^{-1}]" << endl;
        f << "Средняя энергия ионов на питаемом электроде       = " << mean_i_energy_pow << " [eV]" << endl;
        f << "Средняя энергия ионов на заземленном электроде    = " << mean_i_energy_gnd << " [eV]" << endl;
        f << "Поток электронов на питаемом электроде            = " << N_e_abs_pow * WEIGHT / ELECTRODE_AREA / (no_of_cycles * PERIOD) << " [m^{-2} s^{-1}]" << endl;
        f << "Поток электронов на заземленном электроде         = " << N_e_abs_gnd * WEIGHT / ELECTRODE_AREA / (no_of_cycles * PERIOD) << " [m^{-2} s^{-1}]" << endl;

        // рассчитывает усредненное по пространству и времени поглощение энергии электронами и ионами


        auto power_e = accumulate(powere_xt.begin(), powere_xt.end(), 0.0) / static_cast<double>(N_XT * N_G);
        auto power_i = accumulate(poweri_xt.begin(), poweri_xt.end(), 0.0) / static_cast<double>(N_XT * N_G);
        f << line << endl;
        f << "Поглощенная мощность, рассчитанная как <j*E>:" << endl;
        f << "Плотность электронной мощности (средняя)          = " << power_e << " [W m^{-3}]" << endl;
        f << "Плотность мощности ионов (средняя)                = " << power_i << " [W m^{-3}]" << endl;
        f << "Общая плотность мощности (средняя)                = " << power_e + power_i << " [W m^{-3}]" << endl;
        f << line << endl;
        f.close();

    }
}

//------------------------------------------------------------------------------------------//
// main                                                                                     //
// аргументы командной строки:                                                              //
// [1]: количество циклов (0 для инициализации)                                             //
// [2]: "m" включает сбор и сохранение данных                                               //
//------------------------------------------------------------------------------------------//

int main(int argc, char* argv[]) {
    setlocale(LC_ALL, "ru");
    cout << ">> Состояние: Старт программы" << endl;
    cout << "**************************************************************************" << endl;
//   cout << ">> Состояние: Copyright (C) 2021 Z. Donko et al." << endl;
//   cout << ">> Состояние: This program comes with ABSOLUTELY NO WARRANTY" << endl;
//   cout << ">> Состояние: This is free software, you are welcome to use, modify and redistribute it" << endl;
//   cout << ">> Состояние: according to the GNU General Public License, https://www.gnu.org/licenses/" << endl;
//   cout << ">> Состояние: **************************************************************************" << endl;

    if (argc == 1) {
        cout << ">> Состояние: ошибка: нужен аргумент starting_cycle, т.е. 0" << endl;
        return 1;
    }
    else {
        vector<string> argList(argv + 1, argv + argc);
        arg1 = stoi(argList[0]);
        if (argc > 2) {
            if (argList[1] == "m") {
                measurement_mode = true;                            // измерения будут сохранены в файл
            }
            else {
                measurement_mode = false;                           // измерения не снимаются
            }
        }
    }
    if (measurement_mode) {
        cout << ">> Состояние: режим измерения: вкл." << endl;
    }
    else {
        cout << ">> Состояние: режим измерения: выкл." << endl;
    }
    set_electron_cross_sections_ar();
    set_ion_cross_sections_ar();
    calc_total_cross_sections();
    //test_cross_sections(); return 1;

    if (arg1 == 0) {
        ifstream file("picdata.bin", std::ios::binary);
        if (file.good()) {
            file.close();
            cout << ">> Состояние: предупреждение: обнаружены данные из предыдущего расчета." << endl;
            cout << "           Чтобы начать новую симуляцию с самого начала, пожалуйста, удалите все выходные файлы перед запуском программы 0 0" << endl;
            cout << "           Чтобы продолжить существующий расчет, пожалуйста, укажите количество циклов для выполнения, например 100" << endl;
            exit(0);
        }
        no_of_cycles = 1;
        cycle = 1;                                        // иницилизация цикла
        init(N_INIT);                                     // затравочные начальные электроны и ионы
        cout << ">> Состояние: запуск цикла инициализации" << endl;
        Time = 0;
        do_one_cycle();
        cycles_done = 1;
    }
    else {
        no_of_cycles = arg1;                              // выполнение указанного количества циклов
        load_particle_data();                             // чтение предыдущей конфигурации из файла
        if (no_of_cycles % 10 == 1) {
            cout << ">> Состояние: запуск " << no_of_cycles << " цикла" << endl;
        }
        else {
            cout << ">> Состояние: запуск " << no_of_cycles << " циклов" << endl;
        }
        for (cycle = cycles_done + 1;cycle <= cycles_done + no_of_cycles;cycle++) { 
            do_one_cycle(); 
        }
        cycles_done += no_of_cycles;
    }
    datafile.close();
    save_particle_data();
    if (measurement_mode) {
        check_and_save_info();
    }
    if (no_of_cycles%10 == 1)
        cout << ">> Состояние: симуляция " << no_of_cycles << " цикла выполнена." << endl;
    else
        cout << ">> Состояние: симуляция " << no_of_cycles << " циклов выполнена." << endl;
}