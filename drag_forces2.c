#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "drag_forces2.h"

extern Gas1 gas1;

void aerodynamical_gas_drag(struct reb_simulation* const r)
{
    double G = r->G;
    const int N = r->N;
    double time=r->t;
    struct reb_particle* particles = r->particles;
    struct reb_particle com = particles[0]; // center of mass in star initially; 
    
    for(int i=1;i<N;i++)
    {
        struct reb_particle* p = &(particles[i]);
        // posição relativa ao centro de massa 
        double dx = p->x-com.x;    
        double dy = p->y-com.y;
        double dz = p->z-com.z;
        double r = sqrt ( dx*dx + dy*dy + dz*dz );
       // double r = sqrt ( dx*dx + dy*dy );          // é o s do artigo.
        // velocidade relativa ao centro de massa
        double dvx = p->vx-com.vx;  // velocidade da partícula
        double dvy = p->vy-com.vy;
        double dvz = p->vz-com.vz;
        //Size of the plantesimal
        double size=p->r;
        
        // temperatura do gas
        //gas.temperature=pow(r,gas.coef_temp);
        //gas.temperature=gas.gas_temperature_1au*gas.temperature;  // gas.gas_teperatura = 120 K
        // então gas1.T_0 = 120 K
        gas1.T = gas1.T_0 * pow( r, -gas1.beta);  // Eq.(2)  gas1.T_0 é a temperatura a 1AU
        // o gas.beta é o gas.coef_temp
        
        // peso molecular do gás hidrogenio
        gas1.p_molecular=2.3*gas1.massa_molar;
        
        // Velocidade kepleriana
        gas1.velocidade_kepleriana = sqrt(G*com.m/r); // vem da expressão: v = omega x r, omega = sqrt(GM/r³), então: v = omega*r*sin(theta), onde theta = pi/2. Portanto, v = omega*r*1 = omega*r
        // substituindo omega, temos: v = sqrt(GM/r³)*r = sqrt(GM/r) * r/r = sqrt(GM/r) ----> gas1.velocidade_kepleriana
        // fazendo vk = omega x r = -omega*yî + omega*xĵ
        // mas: ê = r/|r| = (x,y,z)/|r| = (x/|r|, y/|r|, z/|r|) = (î, ĵ, k)
        // portanto: vk_x = -omega*y, vk_y = omega*x
        
        // --------------------------------------------------------------- OK!!! -------------------------------------------------------
        
        // O que é gas.H? gas.H é o z_s do artigo
        // Scale Height of the gas: 
        //if(gas.scale_height_1au == 0)  // gas.scale_height_1au é o 0.047 que aparece na Eq.(7)
        //{
        // printf("\n warning - using gas aspect ratio to calculate Scale Height of the gas\n"); 
        // gas.H = gas.gas_aspect_ratio*r; 
        //}
        //else gas.H= gas.scale_height_1au*pow(r,gas.coef_scale_height); // gas.coef_scale_height é o 5/4 que aparece na Eq.(7) ---> Brasser et al. (2007)  adotam esse
        
        // Então:
        
        double cte_H = 0.047;
        double expo_zs = 1.2;
        gas1.z_s = cte_H * pow(r, expo_zs);
        
        // --------------------------------------------------------------- OK!!! -------------------------------------------------------
        
        // Agora vamos definir a velocidade do som. 
        // O artigo de Brasser et al. (2007) usa a definição que está no else:
        
        // Sound speed of gas: 
       //if(gas.sound_speed_1au == 0)
       // {
        // printf("\n warning - using gas.H to calculate sound speed\n"); 
        // gas.sound_speed = gas.velocity_keplerian*r*gas.H; 
         // Implement - second option maybe a flag here: 
         //gas.sound_speed = ( stefan_boltzmann*gas.temperature ) / ( gas.molecular_weight* gas.atom_mass ); 
         //gas.sound_speed = sqrt(gas.sound_speed);
       // }
      //  else
       // {
       // gas.sound_speed=pow(r,gas.coef_sound_speed)*gas.velocity_keplerian; // AU/year  // gas.coef_sound_speed é o 
       // gas.sound_speed=gas.sound_speed_1au*gas.sound_speed; // AU/year
       // }
       
       // A equação para a velocidade do som aparece na Eq.(12) do artigo.
       double coef_sound = -gas1.beta/2.;
       gas1.c_s = 1.43e3 * pow(r, coef_sound);
        
       // --------------------------------------------------------------- OK!!! -------------------------------------------------------
       
       // hora de calcular a densidade 
       // o artigo apresenta duas formas: a primeira que não usa a gaussiana, que é dada pela Eq.(2)
       // a segunda é a que usa a gaussiana, dada pela Eq.(5)
       // usando a Eq.(7)
       
       double dz2 = -dz*dz;
       double z_s2 = gas1.z_s * gas1.z_s;
       gas1.rho_g = gas1.rho0_g * pow(r, gas1.alpha) * exp(dz2/z_s2);  // onde gas1.H é o z_s
       
       // --------------------------------------------------------------- OK!!! -------------------------------------------------------
       
       // eta é dado pela Eq.(4), que é a forma geral, e dado pela Eq.(9), que é o caso específico.
       // fazendo para a Eq.(9):
       gas1.eta = 2.13e-3*pow(r, 0.5);
       double eta1;
       double razao_vel;
       double soma_coef;
       eta1 = 1./(2.*gas1.gamma);
       razao_vel = gas1.c_s/gas1.velocidade_kepleriana;
       soma_coef = gas1.alpha + gas1.beta;
       //gas1.eta   = eta1 * (soma_coef) * pow( razao_vel, 2.);
       
       // --------------------------------------------------------------- OK!!! -------------------------------------------------------
       
       // agora vamos definir a força de arrasto, que é dada pela Eq.(10)
       // mas antes, precisamos escrever a expressão para a velocidade relativa. Ela é dada pela diferença entre a velocidade da partícula e a velocidade do gás.
       // portanto:
       
       double pre1;
       pre1 = 1 - gas1.eta;
       //gas1.v_g   = gas1.v_k * pow(1 - 2 * gas1.eta, 0.5);  // Eq.(3)
       gas1.v_gx  = (-sqrt(pre1)*gas1.velocidade_kepleriana * dy)/r;  // velocidade do gás em x, em y e em z (abaixo)
       gas1.v_gy  = (sqrt(pre1)*gas1.velocidade_kepleriana * dx)/r;
       gas1.v_gz  = 0;
       gas1.v_px  = dvx;  // velocidade da partícula em x, em y e em z (abaixo)
       gas1.v_py  = dvy;
       gas1.v_pz  = dvz;
       
       double v_g;
       
       v_g = gas1.v_gx*gas1.v_gx + gas1.v_gy*gas1.v_gy + gas1.v_gz*gas1.v_gz;
       gas1.v_gas  = sqrt(v_g);
        
        //gas1.eta   = ( 1/(2 * gas1.gamma) ) * (gas1.alpha + gas1.beta) * pow( gas1.c_s/gas1.velocidade_kepleriana, 2);  // Eq.(4)
       // gas1.s     = pow ( (r->x * r->x) + (r->y * r->y), 2 );
       // gas1.z_s   = 0.047 * ( pow( s, 5/4) ); // Eq.(7)
       gas1.v_rx  = gas1.v_px - gas1.v_gx;  // velocidade relativa em x
       gas1.v_ry  = gas1.v_py - gas1.v_gy;  // velocidade relativa em y
       gas1.v_rz  = gas1.v_pz;  // velocidade relativa em z // onde gas.v_kz = 0
        
       gas1.v_r   = gas1.v_rx * gas1.v_rx + gas1.v_ry * gas1.v_ry + gas1.v_rz * gas1.v_rz;   // é o v_r² da equação
       gas1.v_rel = sqrt( gas1.v_r );
       
       // agora precisamos calcular os números adimensionais: Ma, Kn, Re:
       
       gas1.Ma = gas1.v_rel/gas1.c_s;
       gas1.Kn = 1.67e-8/(gas1.rho_g * size);
       gas1.Re = 4.44 * gas1.Ma/gas1.Kn;
       
       double F_D;
       double coef_drag;
       double size2;
       //double density_p;
       size2=size*size;
       coef_drag = drag_coeff(gas1.Ma, gas1.Kn, gas1.Re);
       F_D = 0.5 * coef_drag * M_PI * size2 * gas1.rho_g * gas1.v_rel;
       F_D = F_D/p->m;  // força por unidade de massa -----> aceleração
       //F_D = (3/8) * coef_drag * gas1.rho_g * (gas1.v_rel * gas1.v_rel)/(size * density_p)  // força por unidade de massa
       p->ax += -gas1.v_rx*F_D;   // aceleração em x
       p->ay += -gas1.v_ry*F_D;   // aceleração em y
       p->az += -gas1.v_rz*F_D;   // aceleração em z
       com = reb_get_com_of_pair(com,particles[i]);
       
      // printf("rho = %e vel kepleriana = %e vel relativa = %e v_gasx = %e drag_coef = %e eta = %e vgas = %e\n", gas1.rho_g, gas1.velocidade_kepleriana, gas1.v_rel, gas1.v_gx, dragcoeff, gas1.eta, gas1.v_gas);
       
        }
}

double drag_coeff(double Ma, double Kn, double Re)  // recebe os valores na expressão de F_D
{
   // o artigo mostra duas situações para chegar ao coef de arrasto: Kn < 1 ou Kn > 1
   double dragcoeff;
   if(Kn < 1.)
   {
      if(Ma >= 1.)
      {
         dragcoeff = 2.0;
      }
      else
      {
         if(Re > 1.e3)
         {
            dragcoeff = 0.44 + 1.56 * Ma * Ma;
         }
         else
         {
            dragcoeff = 2. * Ma * Ma + (1. - Ma * Ma) * (2.4e1*(1 + 0.15 * pow(Re, 0.687) )/Re);
         }
      }
   }
   else
   {
      if(Ma < 1.)
      {
         dragcoeff = 3.6/Ma;
      }
      else
      {
         dragcoeff = 2.0;
      }
   }
   
   //printf("rho = %e vel kepleriana = %e vel relativa = %e v_gasx = %e drag_coef = %e eta = %e vgas = %e\n", gas1.rho_g, gas1.velocidade_kepleriana, gas1.v_rel, gas1.v_gx, dragcoeff, gas1.eta, gas1.v_gas);
   return(dragcoeff);
}
