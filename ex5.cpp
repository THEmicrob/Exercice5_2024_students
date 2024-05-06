#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include "ConfigFile.tpp"
#include <algorithm>

using namespace std;

// possibilités pour les conditions aux bords : fixe, libre, sortie, sinusoidale et periodique

void boundary_condition(vector<double> &fnext, vector<double> &fnow, double const &A,
                        double const &t, double const &dt,
                        vector<double> &beta2, string &bc_l, string &bc_r, int &N)
{
  // xL and xR are the left and right positions we have to initialize : indexes 0 and N-1
  if (bc_l == "fixe")
  {
    fnext[0] = fnow[0];
  }
  else if (bc_l == "libre")
  {
    fnext[0] = fnow[1];
    // fnext[0] = fnext[1];
  }
  else if (bc_l == "sortie")
  {
    fnext[0] = fnow[0] + sqrt(beta2[0]) * (fnow[1] - fnow[0]);
  }
  else if (bc_l == "periodique")
  {
    fnext[0] = fnow[N - 1];
  }
  else if (bc_l == "sinusoidale")
  {
    fnext[0] = A * sin(2 * M_PI * t);
  }
  else
  {
    std::cout << "Error: boundary condition not recognized" << endl;
    exit(1);
  }
  if (bc_r == "fixe")
  {
    fnext[N - 1] = fnow[N - 1];
  }
  else if (bc_r == "libre")
  {
    fnext[N - 1] = fnow[N - 2];
    // fnext[N - 1] = fnext[N - 2];
  }
  else if (bc_r == "sortie")
  {
    fnext[N - 1] = fnow[N - 1] - sqrt(beta2[N - 1]) * (fnow[N - 1] - fnow[N - 2]);
  }
  else if (bc_r == "periodique")
  {
    fnext[N - 1] = fnow[0];
  }
  else if (bc_r == "sinusoidale")
  {
    fnext[N - 1] = A * sin(2 * M_PI * t);
  }
  else
  {
    std::cout << "Error: boundary condition not recognized" << endl;
    exit(1);
  }
}

double finit(double x, double xL, double n_init, double xR, double A, double u, string initialization)
{
  double finit_(0.);
  const double PI = 3.1415926535897932384626433832795028841971e0;
  if (initialization == "mode")
  {
    finit_ = 2 * A * sin((n_init + 0.5) * PI * x / (xR - xL)) * cos((n_init + 0.5) * PI * 0 / (xR - xL));
  }
  else if (initialization == "eq4")
  {
    if (x <= xL or x >= xR)
    {
      finit_ = 0.;
    }
    else
    {
      finit_ = 0.5 * A * (1 - cos(2 * PI * (x - xL) / (xR - xL)));
    }
  }
  else
  {
    std::cout << "Error: initialization not recognized" << endl;
    exit(1);
  }
  return finit_;
}

// Surcharge de l'opérateur pour écrire les éléments d'un tableau dans un fichier //

template <class T>
ostream &operator<<(ostream &o, vector<T> const &v)
{
  unsigned int len(v.size());
  for (unsigned int i(0); i < (len - 1); ++i)
    o << v[i] << " ";
  if (len > 0)
    o << v[len - 1];
  return o;
}

// Main //

int main(int argc, char *argv[])
{
  const double PI = 3.1415926535897932384626433832795028841971e0;
  const double g = 9.81;
  double dx;
  double dt;
  double t;
  double Nsteps;
  int stride(0);

  string inputPath("input_example"); // Fichier par défaut
  if (argc > 1)                      // Fichier d'input spécifié par l'utilisateur ("./Exercice7 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockés dans une "map" de strings.

  for (int i(2); i < argc; ++i) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Paramètres de simulation :
  double tfin = configFile.get<double>("tfin");
  int nx = configFile.get<int>("nx"); // nb d'intervalles
  int N = nx + 1;                     // nb de pts de maillage
  double CFL = configFile.get<double>("CFL");
  double nsteps = configFile.get<double>("nsteps");
  double A = configFile.get<double>("A");
  double n_init = configFile.get<double>("n_init"); // numéro du mode propre
  double hL = configFile.get<double>("hL");
  double hR = configFile.get<double>("hR");
  double hC = configFile.get<double>("hC");
  double h00 = configFile.get<double>("h00"); // profondeur, cas uniforme
  double x1 = configFile.get<double>("x1");
  double x2 = configFile.get<double>("x2");
  double xa = configFile.get<double>("xa");
  double xb = configFile.get<double>("xb");
  double xc = configFile.get<double>("xc");
  double xd = configFile.get<double>("xd");
  double xL = configFile.get<double>("xL");
  double xR = configFile.get<double>("xR");
  int n_stride(configFile.get<int>("n_stride"));

  // Conditions aux bords: fixe, libre, sortie de l'onde
  string bc_l = configFile.get<string>("cb_gauche");
  string bc_r = configFile.get<string>("cb_droite");
  // vector<string> bc = {configFile.get<string>("cb_gauche"), configFile.get<string>("cb_droite")};

  // Type de forme initiale de la vague: selon donnée Eq.(4) ou mode propre
  // (par exemple 'mode' pour mode propre, autrement Eq.(4))
  string initialization = configFile.get<string>("initialization");

  // Onde partant initialement vers la gauche ou vers la droite ou statique
  // (par exemple 'gauche', 'droite', 'statique')
  string initial_state = configFile.get<string>("initial_state");

  // Selecteur pour le cas h0 uniforme:
  bool v_uniform = configFile.get<bool>("v_uniform");

  // Sélecteur pour choisir le pas de temps
  bool impose_nsteps = configFile.get<bool>("impose_nsteps");
  // true : dt = tfin/nsteps
  // false : dt = CFL*dx/sqrt(g*h0)

  vector<double> h0(N);   // profondeur aux points de maillage
  vector<double> vel2(N); // u^2 = g*h_0 aux points de maillage
  vector<double> x(N);    // positions des points de maillage
  vector<double> fpast(N), fnow(N), fnext(N), beta2(N);

  dx = (xR - xL) / (N - 1);
  bool ecrire_f = configFile.get<bool>("ecrire_f"); // Exporter f(x,t) ou non

  // Eq.(1) ou Eq.(2) [ou Eq.(6) (facultatif)]: Eq1, Eq2 ou Eq6
  string equation_type = configFile.get<string>("equation_type");

  // Initialisation de la profondeur et de la vitesse aux points de maillage : //
  for (int i(0); i < N; ++i)
  {
    // Initialisation de la position //
    x[i] = xL + i * dx;
    // Initialisation de la profondeur //
    if (v_uniform)
    {
      h0[i] = h00;
    }
    else
    { // profondeur dépend de la position
      if (x[i] <= xa and x[i] >= xL)
      {
        h0[i] = hL;
      }
      else if (x[i] < xb)
      {
        h0[i] = 0.5 * (hL + hC) + 0.5 * (hL - hC) * cos(PI * (x[i] - xa) / (xb - xa));
      }
      else if (x[i] <= xc)
      {
        h0[i] = hC;
      }
      else if (x[i] < xd)
      {
        h0[i] = 0.5 * (hC + hR) - 0.5 * (hR - hC) * cos(PI * (x[i] - xc) / (xd - xc));
      }
      else if (x[i] <= xR)
      {
        h0[i] = hR;
      }
      else
      {
        std::cout << "Error: x[i] out of bounds" << endl;
        exit(1);
      }
    }
    // Initialisation de la vitesse au carré //
    vel2[i] = g * h0[i];
  }
  // TO CHECK

  auto max_vel2 = std::max_element(vel2.begin(), vel2.end());

  // Initialisation du pas de temps //
  dt = CFL * dx / sqrt(*max_vel2);
  if (impose_nsteps)
  { // redéfinition de dt pour que CFL=1
    dt = tfin / nsteps;
    CFL = dt * sqrt(*max_vel2) / dx;
  }
  // TO CHECK

  // Fichiers de sortie :
  string output = configFile.get<string>("output");

  ofstream fichier_x((output + "_x").c_str());
  fichier_x.precision(15);

  ofstream fichier_v((output + "_v").c_str());
  fichier_v.precision(15);

  ofstream fichier_f((output + "_f").c_str());
  fichier_f.precision(15);

  // Initialisation des tableaux du schéma numérique : //
  for (int i(0); i < N; ++i)
  {
    // Initialisation de fnow //
    fnow[i] = finit(x[i], x1, n_init, x2, A, sqrt(vel2[i]), initialization);

    // Initialisation de fpast //
    if (initial_state == "statique")
    { // système initialement au repos
      fpast[i] = fnow[i];
    }
    else if (initial_state == "gauche")
    { // propagation rétrograde (vers la gauche)
      fpast[i] = finit(x[i] - sqrt(vel2[i]) * dt, x1, n_init, x2, A, sqrt(vel2[i]), initialization);
    }
    else if (initial_state == "droite")
    { // propagation progressive (vers la droite)
      fpast[i] = finit(x[i] + sqrt(vel2[i]) * dt, x1, n_init, x2, A, sqrt(vel2[i]), initialization);
    }
    else
    {
      std::cout << "Error: initial_state not recognized" << endl;
      exit(1);
    }

    // Initialisation de beta^2 //
    beta2[i] = vel2[i] * dt * dt / (dx * dx);
  }
  // TO CHECK

  cout << "beta2[0] is " << beta2[0] << endl;
  cout << "dt is " << dt << endl;

  // Boucle temporelle :
  for (t = 0.; t < tfin - .5 * dt; t += dt)
  {
    // Ecriture :
    if (stride % n_stride == 0)
    {
      if (ecrire_f)
        fichier_f << t << " " << fnow << endl;
    }
    ++stride;

    // Evolution :
    for (int i(1); i < N - 1; ++i)
    {
      if (equation_type == "Eq1")
      {
        // Eq.(1) D2F/DT2 = D/DX(g*h0*DF/DX)
        fnext[i] = 2 * (1 - beta2[i]) * fnow[i] - fpast[i] + beta2[i] * (fnow[i + 1] + fnow[i - 1]) * (1 + (h0[i + 1] - h0[i - 1]) / (2 * h0[i]));
        // fnext[i] = 2*(1-beta2[i])*fnow[i] - fpast[i] + beta2[i]*(fnow[i+1]+fnow[i-1])+beta2[i]/h0[i]*(h0[i+1]-h0[i])*(fnow[i+1]-fnow[i]);
        // TO CHECK
      }
      else if (equation_type == "Eq2")
      {
        // Eq.(2) D2F/DT2 = g*h0*D2F/DX2
        fnext[i] = 2 * (1 - beta2[i]) * fnow[i] - fpast[i] + beta2[i] * (fnow[i + 1] + fnow[i - 1]);
      }
      else if (equation_type == "Eq6")
      {
        // Eq.(6) D2F/DT2 = D2/DX2(g*h0*F)
        fnext[i] = 2 * (1 - beta2[i]) * fnow[i] - fpast[i] + beta2[i] * (fnow[i + 1] + fnow[i - 1]);
        // totally wrong
      }
      else
      {
        std::cout << "Error: equation_type not recognized" << endl;
        exit(1);
      }
    }

    // Conditions aux bords //
    boundary_condition(fnext, fnow, A, t, dt, beta2, bc_l, bc_r, N);

    // Mise à jour des tableaux //
    fpast = fnow;
    fnow = fnext;
  }

  if (ecrire_f)
  fichier_f << t << " " << fnow << endl;
  fichier_x << x << endl;
  fichier_v << vel2 << endl;

  fichier_f.close();
  fichier_x.close();
  fichier_v.close();

  return 0;
}
