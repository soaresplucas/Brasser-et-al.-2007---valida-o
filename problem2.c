#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "tools.h"
//#include "../../forces/drag_forces.h"
#include "../../forces/drag_forces2.h"
#include <string.h>
#define LENGHT 400


const double G =  6.67435e-11; // cte universal da gravitção
const double ua = 149597870700;  // unidade astronomica (em metros)
const double grams = 1e-3; //kg
const double cm = 6.68458712e-14; //au - 1/(au*100)
const double M_sun=1.98841583e30; // kg - Resolution IAU
const double M_earth=5.9722e24; // kg
const double year=31557600; // s
const double km=1e3; //m
const double m_cm=100.; //cm
const double avogrado_number=6.02214076e23; // g/mol
const double sun_radius=0.00465047; // AU - (6.957 x 10^8m / 1 au in m) 
const int N_planetesimais=2;
const double segundos = 3.17098e-8; // 1 segundo = 3.17098e-8 anos


double ecc[2]=
{
    0.,
    0.2,
    //0.,
};

double I[2]= //degrees 
{
    10.,
    0.,
};


Gas1 gas1;

void heartbeat(struct reb_simulation* const r);


// multipurpose structure
typedef struct Normalized 
{
 double N_elements;
 double tmax;
} Normalized;
Normalized normal;


int main(int argc, char* argv[])
{
    double cm3=cm*cm*cm;
    struct reb_simulation* r = reb_create_simulation();
    r->integrator= REB_INTEGRATOR_IAS15;
    r->force_is_velocity_dependent = 1; // Force only depends on velocities.
   // r->ri_ias15.epsilon             = 1e-4;
    r->G=G/ua/ua/ua;
    r->G=r->G*year*year;
    r->dt=1e-6;
    struct reb_particle sun = {0};
    sun.m=M_sun;
    sun.r=sun_radius;
    sun.hash = reb_hash("Sun");
    reb_add(r, sun);
    struct reb_particle primary = {0}; 
    primary = r->particles[0];
    
    //Gas properties - structure Hydrogeny 
    
    gas1.T_0 = 280.;  //
    gas1.beta = 0.5; //
    gas1.massa_molar = 2.01588;  //
    
    // Eq.(2)
    //gas1.rho0_g = 1.4e-9;  // ver as unidades
    //gas1.rho0_g = gas1.rho0_g/cm3;
    //gas1.rho0_g = gas1.rho0_g*grams;

    gas1.rho0_g = (1.4e-9)*(grams/cm3);
    //gas1.s      = 
    gas1.alpha  = -11./4.;
    //gas1.T_0     = 280;  // K
    //gas1.beta   = -1/2;
    
    //gas1.massa_molar = 2.01588;
    gas1.massa_molar = gas1.massa_molar*grams;
    
    gas1.gamma = 7./5.;
    
    //gas1.R = 8.314;  // J K-1 mol-1
    //gas1.R = 8.314e3;  // g (m/s)² K-1 mol-1
    gas1.R = 8.314e3 * (100. * cm)*(100. * cm)/(segundos*segundos);  // g ( AU/yr ) K-1 mol-1
    
    //gas1.rho_g = gas1.rho_g/cm3;
    
    r->N_active = r->N;
    

    double m;
    double a,exc,inc,Omega,omega,f;
    double size,density_p;
    double cm_m=1e-2;
    size=2*km;
    density_p=2.0; // g/cm^3
    density_p=density_p/(cm_m*cm_m*cm_m); // g/m^3
    m=density_p*(4./3.)*M_PI*size*size*size; // g
    m=m*grams; //kg
    a=1.;
    f=0.;
    omega=0.;
    Omega=0.;
    size=size/ua; //  UA conversion
    r->dt=1e-6;
    
    
    for(int j=0;j<N_planetesimais;j++)
    {
        exc=ecc[j];
        inc=I[j]*M_PI/180.; // rad 
        struct reb_particle planetesimais = reb_tools_orbit_to_particle(r->G, primary, m, a, exc, inc, Omega, omega, f);
        planetesimais.r=size;
        planetesimais.hash = j+1;
        reb_add(r, planetesimais);
        //printf("%f\n", F_D);
    }
    normal.tmax=2e4;
    reb_move_to_com(r);   
    r->additional_forces= aerodynamical_gas_drag;  
    r->heartbeat = heartbeat;
    // Start integration
    reb_integrate(r, normal.tmax);   // começa a integração
    printf("Done.\n");
    
    printf("Re = %f\n", gas1.Re);
    printf("acabou!!!!!!!!!!!!!!!!!!!!!\n");
    
}    

void heartbeat(struct reb_simulation* const r)
{
  if (reb_output_check(r, 1000)) reb_output_timing(r, normal.tmax);
    int i;
    double intervaltime2=10;
    int timecontrol2=reb_output_check(r,intervaltime2);
    char namefile2[100]; 
    int erro=0;
    struct reb_particle com = r->particles[0];
    double tempoanos=r->t;// planets: 
    for(i=1;i<r->N;i++)
    {
        struct reb_orbit part = reb_tools_particle_to_orbit_err(r->G, r->particles[i],com,&erro);
        if (timecontrol2)
        {
            snprintf (namefile2,sizeof(namefile2), "part-%.d.out", i);
            FILE*arq2;
            arq2 = fopen(namefile2,"a");
            fprintf(arq2,"%e %e %e %e  \n",tempoanos,part.a,part.e,part.inc);
            fclose(arq2);
        }
    }
}    
    
