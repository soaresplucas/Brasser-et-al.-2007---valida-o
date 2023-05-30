/**
 * Brasser et al. (2007):  
 * Erikson et al. (2021); Perets and Murray-Clay, 2011; Guillot et al. (2014): 
 * Dynamical friction force (Grishin and Perets (2015)): 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include "rebound.h"

// Classe Gas
typedef struct Gas1
{
  // ATRIBUTOS
  
  double p_molecular;
  double massa_molar;
  
  // Eq.(1)
  double v_gx; // velocidade do gás em x 
  double v_gy; // velocidade do gás em y
  double v_gz; // velocidade do gás em z
  double v_kx; // velocidade kepleriana em x do gás
  double v_ky; // velocidade kepleriana em y do gás
  double v_kz; // velocidade kepleriana em z do gás
  double velocidade_kepleriana;
  double w;  // velocidade angular do gás (para o calculo da velocidade kepleriana)
  double rho_g; // densidade do gás
  double P_g; // pressão do gás
  double H; // scale height
  double v[3]; // velocities
  double gaslifetime;
  double v_gas;
  
  // Eq.(2)
  double rho0_g; // densidade do gás a 1 AU
  double T_0; // temperatura do gás a 1 AU
  double T; // temperatura do gás
  
  // Eq.(3)
  // v_g = v_k * (1 - 2*eta)^1/2
  
  // Eq.(4)
  double eta;
  double c_s; // velocidade do som
  double gamma; // razão entre as capacidades termicas especificas do gás  --> dado após a Eq.(18): gamma = 7/5
  double alpha;
  double beta;
  
  // Eq.(6)
  double s; // s² = x² + y² 
  
  // Eq.(7)
  double z_s;
  
  // Eq.(10)  -->  equação para a força de arrasto
  double C_D; // coeficiente de arrasto
  double Re;  // numero de Reynolds
  double r_c; // raio do cometa
  double v_rx; // velocidade relativa em x entre o cometa e o gás
  double v_ry; // velocidade relativa em y entre o cometa e o gás
  double v_rz; // velocidade relativa em z entre o cometa e o gás
  double v_r;   // componentes ao quadrado da velocidade relativa
  double v_rel; // velocidade relativa
  double v_px; // velocidade da partícula em x para o calculo da velocidade relativa v_r
  double v_py; // velocidade da partícula em y para o calculo da velocidade relativa v_r
  double v_pz; // velocidade da partícula em z para o calculo da velocidade relativa v_r
  
  // subseção 2.1
  double Ma; // numero de Mach
  double Kn; // numero de knudsen
  double lambda;  // livre caminho médio
  
  // Eq.(12)
  double R; // cte dos gases
  double m_H; // massa de um mol de hidrogenio
  double e; // excentricidade -> no rebound tem como obter este valor
  
  // Eq.(13)  --> não tem novos atributos
  // Eq.(14)  --> não tem novos atributos
  // Eq.(15)
  double visc;  // viscosidade do gás
  
  // Eq.(16)
  double c_m; // velocidade molecular
  
  // Eq.(17)  --> calculo de c_m
  // Eq.(18)  --> não tem novos atributos
  
  // Eq.(19)
  double m_p; // massa do proton
  double CS; // corte transversal da molecula de hidrogenio CS = pi*d²
  double d;  // diametro da meloceula de hidrogenio
  
  // Eq.(20)  --> não tem novos atributos
  
  // Eq.(21)  --> é a eq para z_c, como estamos no plano, não precisa levar em conta
  double s_c;
  
  // Eq.(22)  --> calcula o coeficiente de arrasto
  // Eq.(23)  -->  calcula a força de arrasto levando em conta a Eq.(22)
  
  int ative_agd; // ative Aerodynamical drag force put > 0 
  int ative_gdf; // ative Gas Friction force > 0   
  int ative_interpol;
  
  
 // double massa_total;
 // double densidade_superficial;
 // double densidade_superficial0;  // a 1 AU
 // double densidade0;  // densidade do gás a 1 AU
 // double densidade;   // densidade do gás
 // double temperatura0;  // tenperatura do gás a 1 AU
 // double temperatura;  // temperatura do gas
 // double peso_molecular;
 // double grad_presssao; // gradiente de pressão
 // double pressao_gas;  // pressão do gás a 1 AU
 // double v_kepleriana; // velocidade kepleriana
 // double v_gas; // velocidade do gas
 // double v_relativa; // velocidade relativa
 // double v_som;
 // double livre_caminho_medio;
 // double massa_molecular;
 // double gamma;  // é o gamma do artigo Brasser et al.(2007)
 // double drag_coefficiente;  // para o calculo da força de arrasto
 // double r_c; // raio da particula
 // double mach_number; // numero de Mach
 // double knudsen_number; // número de Knudsen
 // double reynolds_number;  // numero de reynolds
 // double viscosidade_molecular;
 // double massa_molar;
 // double gaslifetime;
 // double coef_temp;
 // double coef_v_som; 
 // double numero_avogrado;
 // double H; // escala de altura
 // double H_0;  // aspect ratio of the gas
  
} Gas1;

extern Gas1 gas1;
void aerodynamical_gas_drag(struct reb_simulation* const r);
double drag_coeff(double mach_number, double knudsen_number,double reynolds_number);
