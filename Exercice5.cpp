#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include "ConfigFile.tpp"
#include <algorithm>

using namespace std;




void boundary_condition(vector<double> &fnext, vector<double> &fnow, double const& A, \
		double const& t,double const& dt, \
		vector<double> &beta2, string &bc_l, string &bc_r, int &N)
{
  // TODO: Insert boundary conditions
  if (bc_l == "fixe") {
	  fnext[0] = fnow[0];
  }
  else if (bc_l == "libre") {
	  fnext[0] = fnext[1];
  }
  else if (bc_l == "sortie") {
	  fnext[0] = fnow[0] - sqrt(beta2[0])*(fnow[0] - fnow[1]);
  }
  
  if (bc_r == "fixe") {
	  fnext[N] = fnow[N];
  }
  else if (bc_r == "libre") {
	  fnext[N] = fnext[N-1];
  }
  else if (bc_r == "sortie") {
	  fnext[N] = fnow[N] - sqrt(beta2[N])*(fnow[N] - fnow[N-1]);
  }
}


double finit(double x, double xL, double n_init, double xR, double A, double u, string init)
{
  double finit_(0.);
  const double PI = 3.1415926535897932384626433832795028841971e0;
  // TOTO initialiser un mode propre
  double t0=0;
  if (init == "mode") {
	  finit_ = 2*A*sin((n_init+0.5)*PI*x/(xR-xL))*cos((n_init+0.5)*PI*t0*u/(xR-xL));
  }
  else if (init == "autre") {
	  if (x <= xL) {
		  finit_ = 0;		
	  }
	  else if (x >= xR) {
		  finit_ = 0;
	  }
	  else {
		  finit_ = 0.5*A*(1 - cos(2*PI*(x-xL)/(xR-xL)));
	  }
	  
  }
  return finit_;
}

//
// Surcharge de l'operateur pour ecrire les elements d'un tableau
//
template <class T> ostream& operator<< (ostream& o, vector<T> const& v)
{
  unsigned int len(v.size());
  for(unsigned int i(0); i < (len - 1); ++i)
    o << v[i] << " ";
  if(len > 0)
    o << v[len-1];
  return o;
}

//
// Main
//
int main(int argc, char* argv[])
{
  const double PI = 3.1415926535897932384626433832795028841971e0;
  const double g  = 9.81;
  double dx;
  double dt;
  double t;
  double Nsteps;
  int stride(0);

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Parametres de simulation :
  double tfin    = configFile.get<double>("tfin");
  int nx          = configFile.get<int>("nx"); // nb d'intervalles
  int N = nx+1;                                // nb de pts de maillage
  double CFL     = configFile.get<double>("CFL");
  double nsteps  = configFile.get<double>("nsteps");
  double A       = configFile.get<double>("A");
  double n_init   = configFile.get<double>("n_init");
  double hL      = configFile.get<double>("hL");
  double hR      = configFile.get<double>("hR");
  double hC      = configFile.get<double>("hC");
  double h00     = configFile.get<double>("h00"); // profondeur, cas uniforme
  double x1      = configFile.get<double>("x1");
  double x2      = configFile.get<double>("x2");
  double xa      = configFile.get<double>("xa");
  double xb      = configFile.get<double>("xb");
  double xc      = configFile.get<double>("xc");
  double xd      = configFile.get<double>("xd");
  double xL      = configFile.get<double>("xL");
  double xR      = configFile.get<double>("xR");
  int n_stride(configFile.get<int>("n_stride"));
 
// Conditions aux bords:
  string bc_l           = configFile.get<string>("cb_gauche");
  string bc_r           = configFile.get<string>("cb_droite");

// Type de forme initiale de la vague: selon donnée Eq.(4) ou mode propre
// (par exemple 'mode' pour mode propre, autrement Eq.(4))
  string initialization = configFile.get<string>("initialization"); 

// Onde partant initialement vers la gauche ou vers la droite ou statique
// (par exemple 'gauche', 'droite', 'statique')
  string initial_state = configFile.get<string>("initial_state");

// Selecteur pour le cas h0 uniforme:
  bool v_uniform        = configFile.get<bool>("v_uniform");

// Selecteur pour choisir le pas de temps:
// true --> dt=tfin/nsteps; t final est exactement tfin
// false --> dt tel que beta_CFL=1; attention, t final n'est pas exactement tfin
  bool impose_nsteps    = configFile.get<bool>("impose_nsteps");

  vector<double> h0(N) ;   // profondeur aux points de maillage
  vector<double> vel2(N) ; // u^2 = g*h_0 aux points de maillage
  vector<double> x(N) ;    // positions des points de maillage
  vector<double> fpast(N), fnow(N), fnext(N), beta2(N);

  dx = (xR - xL) / (N-1);
  bool ecrire_f = configFile.get<bool>("ecrire_f"); // Exporter f(x,t) ou non

 // Eq.(1) ou Eq.(2) [ou Eq.(6) (facultatif)]: Eq1, Eq2 ou Eq6
  string equation_type = configFile.get<string>("equation_type");  

  for(int i(0); i<N; ++i){ 
    // TODO initialize the depth h0 and the velocity^2 vel2
    x[i] = xL + i*dx;
    
  if (v_uniform){
		h0[i] = h00;
	}
	else {
      if ((xL <= x[i]) & (x[i] <= xa)) {
        h0[i] = hL;
      }
      else if ((xa < x[i]) & (x[i] < xb)) {
        h0[i] = 0.5*(hL+hC) + 0.5*(hL-hC)*cos(PI*(x[i]-xa)/(xb-xa));
      }
      else if ((xb <= x[i]) & (x[i] <= xc)) {
        h0[i] = hC;
      }
      else if ((xc < x[i]) & (x[i] < xd)) {
        h0[i] = 0.5*(hR+hC) + 0.5*(hR-hC)*cos(PI*(x[i]-xd)/(xd-xc));
      }
      else if ((xd <= x[i]) & (x[i] <= xR)) {
        h0[i] = hR;
      }
    }
	
	vel2[i] = g*h0[i];
  }

  auto max_vel2 = std::max_element(vel2.begin(), vel2.end());

  // TODO
  // define the dt according to CLF input 
  if(impose_nsteps){
    // define the dt and CLF when you want to fix nsteps
    dt  = tfin/nsteps;
    // CFL = double(*max_vel2)*dt/dx;
  }
  else {
	  dt = CFL*dx/double(*max_vel2);
	  //tfin = nsteps*dt;
  }

  // Fichiers de sortie :
  string output = configFile.get<string>("output");

  ofstream fichier_x((output.substr(0, output.size()-4) + "_x.out").c_str());
  fichier_x.precision(15);

  ofstream fichier_v((output.substr(0, output.size()-4) + "_v.out").c_str());
  fichier_v.precision(15);

  ofstream fichier_f((output.substr(0, output.size()-4) + "_f.out").c_str());
  fichier_f.precision(15);


  // Initialisation des tableaux du schema numerique :

  //TODO initialize f and beta^2
  for(int i(0); i<N; ++i)
  {
	if (initial_state == "gauche") {
		fpast[i] = finit(x[i] - g*h0[i]*dt, x1, n_init, x2, A, g*h0[i], initialization);
	}
	else if (initial_state == "droite") {
		fpast[i] = finit(x[i] + g*h0[i]*dt, x1, n_init, x2, A, g*h0[i], initialization);
	}
	else if (initial_state == "static") {
		fpast[i] = finit(x[i], x1, n_init, x2, A, g*h0[i], initialization);
	}
    fnow[i]  = finit(x[i], x1, n_init, x2, A, g*h0[i], initialization);
    beta2[i] = pow(vel2[i]*dt/dx, 2);

    // TODO initialize beta2, fnow and fpast according to the requests
  }


  cout<<"beta2[0] is "<<beta2[0]<<endl;
  cout<<"dt is "<< dt <<endl;


  // Boucle temporelle :
  for(t=0.; t<tfin-.5*dt; t+=dt)
  {
    // Ecriture :
    if(stride%n_stride == 0)
    {
      if(ecrire_f) fichier_f << t << " " << fnow << endl;
    }
    ++stride;

    // Evolution :
    for(int i(1); i<N-1; ++i)
    {
      // TODO: write the expressions for fnext 
      if (equation_type=="Eq1")
      {
        fnext[i] = 2*(1-beta2[i])*fnow[i] - fpast[i] + beta2[i]*(fnow[i+1]+fnow[i-1])+beta2[i]/h0[i]*(h0[i+1]-h0[i])*(fnow[i+1]-fnow[i]);
      }
      else if (equation_type=="Eq2")
      {
        fnext[i] = 2*(1-beta2[i])*fnow[i] - fpast[i] + beta2[i]*(fnow[i+1]+fnow[i-1]);
      }
      
    }

    // TODO add boundary conditions
    boundary_condition(fnext, fnow, A, t, dt, beta2, bc_l, bc_r, N);

    //TODO: faire la mise a jour :
    fpast = fnow ;
    fnow  = fnext ;
  }

  if(ecrire_f) fichier_f << t << " " << fnow << endl;
  fichier_x << x << endl;
  fichier_v << vel2 << endl;

  fichier_f.close();
  fichier_x.close();
  fichier_v.close();

  return 0;
}
