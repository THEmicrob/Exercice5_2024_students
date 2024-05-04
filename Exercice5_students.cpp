#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include "ConfigFile.tpp"
#include <algorithm>

using namespace std;

// 3 possibilités pour les conditions aux bords : fixe, libre et sortie de l'onde //

void boundary_condition(vector<double> &fnext, vector<double> &fnow, double const &A,
                        double const &t, double const &dt,
                        vector<double> &beta2, vector<string> &bc, int &N)
{
  for (int i(0); i < 1; ++i)
  {
    double x = 0.0;
    if (i > 0)
    {
      x = N - 1;
    }
    // xL and xR are the left and right positions we have to initialize : indexes 0 and N-1

    if (bc[i] == "fixe")
    {
      fnext[x] = 0.0;
    }
    else if (bc[i] == "libre")
    {
      fnext[x] = fnow[x + 1];
    }
    else if (bc[i] == "sortie")
    {
      if (i == 0)
      {
        fnext[x] = fnow[x] + sqrt(beta2[x]) * (fnow[x + 1] - fnow[x]);
      }
      else
      {
        fnext[x] = fnow[x] - sqrt(beta2[x]) * (fnow[x] - fnow[x - 1]);
      }
    }
    else if (bc[i] == "periodique")
    {
      if (i == 0)
      {
        fnext[x] = fnow[N - 1];
      }
      else
      {
        fnext[x] = fnow[0];
      }
    }
    else if (bc[i] == "sinusoidale")
    {
      fnext[x] = A * sin(2 * M_PI * t);
    }
    else
    {
      std::cout << "Error: boundary condition not recognized" << endl;
      exit(1);
    }
  }
}
// next week

double finit(double x, double x1, double n_init, double x2, string initialization, double A = 1.)
{
  double finit_(0.);
  const double PI = 3.1415926535897932384626433832795028841971e0;
  if (initialization == "mode")
  {
    // TO CHECK
    // Check the constant C : comes from boundary conditions
    double C = 1.;
    finit_ = C * cos(2*PI*n_init*(x-x1)/(x2-x1));
  }
  else if (initialization == "Eq4")
  {
    if (x < x1 or x > x2)
    {
      std::cout << "Error: x out of bounds" << endl;
      exit(1);
    }
    finit_ = 0.5 * A * (1 - cos(2 * PI * (x - x1) / (x2 - x1)));
  }
  else
  {
    std::cout << "Error: initialization not recognized" << endl;
    exit(1);
  }

  return finit_;
}
// next week : partie analytique me manque

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
  // string bc_l = configFile.get<string>("cb_gauche");
  // string bc_r = configFile.get<string>("cb_droite");
  vector<string> bc = {configFile.get<string>("cb_gauche"), configFile.get<string>("cb_droite")};

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
    fnow[i] = finit(x[i], x1, n_init, x2, initialization, A);

    // Initialisation de fpast //
    if (initial_state == "statique")
    { // système initialement au repos
      fpast[i] = fnow[i];
    }
    else if (initial_state == "gauche")
    { // propagation rétrograde (vers la gauche)
      fpast[i] = finit(x[i] - (*max_vel2) * dt, x1, n_init, x2, initialization, A);
    }
    else if (initial_state == "droite")
    { // propagation progressive (vers la droite)
      fpast[i] = finit(x[i] + (*max_vel2) * dt, x1, n_init, x2, initialization, A);
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
      double DF = (fnow[i + 1] - fnow[i - 1]) / (2 * dx);
      double D2F = (fnow[i + 1] - 2 * fnow[i] + fnow[i - 1]) / (dx * dx);
      // fnext[i] = 0.0;
      if (equation_type == "Eq1")
      {
        // Eq.(1) D2F/DT2 = D/DX(g*h0*DF/DX)
      }
      else if (equation_type == "Eq2")
      {
        // Eq.(2) D2F/DT2 = g*h0*D2F/DX2
      }
      else if (equation_type == "Eq6")
      {
        // Eq.(6) D2F/DT2 = D2/DX2(g*h0*F)
      }
      else
      {
        std::cout << "Error: equation_type not recognized" << endl;
        exit(1);
      }
    }

    // Conditions aux bords //
    boundary_condition(fnext, fnow, A, t, dt, beta2, bc, N);

    // Mise à jour des tableaux //
    for (int i(0); i < N; ++i)
    {
      fpast[i] = fnow[i];
      fnow[i] = fnext[i];
    }
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
