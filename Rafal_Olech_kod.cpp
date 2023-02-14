#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

double detJTab[4] = { 0,0,0,0 };


struct node
{
    double x; // współrzędna węzła na osi x
    double y; // współrzędna węzła na osi y
    double temp0 = 100; // temperatura poczatkowa w wezle 
};

struct GlobalData
{
    int i, j, nPC = 4; // liczba punktów całkowania
    int t0 = 1200;     // ambient temperature [C] - temperatura otoczenia 
    int initialTemp = 100; // initial temperature - temperatura początkowa
    int stime = 500;   // simulation time - czas symulacji
    double sst = 50.;  // simulation step time - czas kroku symulacji
    //double sst = 1.;  // simulation step time - czas kroku symulacji - dla testów siatki 31x31
    double c = 700.;   // ciepło właściwe
    double ro = 7800.; // gęstość
    double conductivity = 25; // przewodność
    double alfa = 300.; // alfa
};

//stuktura opisująca punkt całkowania
struct point
{
    double ksi;
    double eta;
    point() {};
    point(double ksi, double eta) : eta(eta), ksi(ksi) {};
};

//kazdy element posiada 4 ściany.
struct wall {
    point punkt1;
    point punkt2;
    bool BC;
    wall() { BC = false; }; // flaga określająca występowanie warunku brzegowego na ścianie domyslnie ustawione jest na false  
    wall(point punkt1, point punkt2)
    {
        this->punkt1 = punkt1;
        this->punkt2 = punkt2;
    };
};


struct element {
    int ID[4];
    wall TabSciany[4];
};

// tablice dynamiczne do pochodnych ksi, eta
struct element4
{
    double** ksi;
    double** eta;

    element4()
    {
        ksi = new double* [4];
        eta = new double* [4];

        for (int i = 0; i < 4; i++)
        {
            ksi[i] = new double[4];
            eta[i] = new double[4];

        }
    }
};


struct grid {
    double H, B;
    int nH, nB, nN, nE;
    node* tabNodes;
    element* tabElements;


    grid(double h = 0.1, double b = 0.1, int nh = 4, int nb = 4)
        : H(h), B(b), nH(nh), nB(nb)
    {
        nE = (nH - 1) * (nB - 1); // liczba elementów
        nN = nH * nB; // liczba węzłów
        tabElements = new element[nE]; // tablica elementów
        tabNodes = new node[nN]; // tablica węzłów


        // fragment kodu odpowiedzialny za obliczanie ID wierzchołków poszczególnych elementów
        int counter = 1; // licznik wspomagający przy obliczaniu ID węzłów 
        for (int i = 0; i < nE; i++)
        {
            element e;
            if (i != 0 && (i + counter) % nH == 0) // warunek obsługujący newralgiczne miejsca numerowania wierzchołków
                counter++;

            e.ID[0] = i + counter;
            e.ID[1] = e.ID[0] + nH;
            e.ID[2] = e.ID[1] + 1;
            e.ID[3] = e.ID[0] + 1;
            tabElements[i] = e;
        }

        // fragment kodu odpowiedzialny za obliczanie współrzędnych (x,y) wierzchołków
        int counterX = 0;
        int counterY = 0;
        for (int i = 0; i < nN; i++) // przejście po węzłach siatki 
        {
            // najpierw wyznaczane są współrzedne wierzchołków po osi Y, 
            // jeżeli w pionie skończą się wierzchołki przeskakuje o 1 w prawo
            if (i != 0 && i % nH == 0)
            {
                counterX++;
            }
            tabNodes[i].x = B / (nB - 1) * counterX; // obliczanie współrzędnej x
            tabNodes[i].y = H / (nH - 1) * counterY; // obliczanie współrzędnej y
            counterY++;
            if (counterY % nH == 0)
            {
                counterY = 0;
            }
        }
    }
};


// struktura z możliwymi współrzędnymi punktów całkowania 
struct pcKsiEta
{
    double a = 1;
    double b = 1 / sqrt(3);
};



struct UniversalElement {

    // współrzędne punktów całkowania ksi, eta dla elementu uniwersalnego (przydaje się do macierzy C)
    double wspolrzednePC1KsiEta[2] = { -(1 / sqrt(3)), -(1 / sqrt(3)) };
    double wspolrzednePC2KsiEta[2] = { (1 / sqrt(3)), -(1 / sqrt(3)) };
    double wspolrzednePC3KsiEta[2] = { (1 / sqrt(3)), (1 / sqrt(3)) };
    double wspolrzednePC4KsiEta[2] = { -(1 / sqrt(3)), (1 / sqrt(3)) };


    double PC1NToTabC[1][4];
    double PC2NToTabC[1][4];
    double PC3NToTabC[1][4];
    double PC4NToTabC[1][4];


    UniversalElement(double x)
    {
        // wyliczanie funkcji kształtu dla każdego punktu całkowania
        PC1NToTabC[0][0] = 0.25 * (1 - wspolrzednePC1KsiEta[0]) * (1 - wspolrzednePC1KsiEta[1]);
        PC1NToTabC[0][1] = 0.25 * (1 + wspolrzednePC1KsiEta[0]) * (1 - wspolrzednePC1KsiEta[1]);
        PC1NToTabC[0][2] = 0.25 * (1 + wspolrzednePC1KsiEta[0]) * (1 + wspolrzednePC1KsiEta[1]);
        PC1NToTabC[0][3] = 0.25 * (1 - wspolrzednePC1KsiEta[0]) * (1 + wspolrzednePC1KsiEta[1]);

        PC2NToTabC[0][0] = 0.25 * (1 - wspolrzednePC2KsiEta[0]) * (1 - wspolrzednePC2KsiEta[1]);
        PC2NToTabC[0][1] = 0.25 * (1 + wspolrzednePC2KsiEta[0]) * (1 - wspolrzednePC2KsiEta[1]);
        PC2NToTabC[0][2] = 0.25 * (1 + wspolrzednePC2KsiEta[0]) * (1 + wspolrzednePC2KsiEta[1]);
        PC2NToTabC[0][3] = 0.25 * (1 - wspolrzednePC2KsiEta[0]) * (1 + wspolrzednePC2KsiEta[1]);

        PC3NToTabC[0][0] = 0.25 * (1 - wspolrzednePC3KsiEta[0]) * (1 - wspolrzednePC3KsiEta[1]);
        PC3NToTabC[0][1] = 0.25 * (1 + wspolrzednePC3KsiEta[0]) * (1 - wspolrzednePC3KsiEta[1]);
        PC3NToTabC[0][2] = 0.25 * (1 + wspolrzednePC3KsiEta[0]) * (1 + wspolrzednePC3KsiEta[1]);
        PC3NToTabC[0][3] = 0.25 * (1 - wspolrzednePC3KsiEta[0]) * (1 + wspolrzednePC3KsiEta[1]);

        PC4NToTabC[0][0] = 0.25 * (1 - wspolrzednePC4KsiEta[0]) * (1 - wspolrzednePC4KsiEta[1]);
        PC4NToTabC[0][1] = 0.25 * (1 + wspolrzednePC4KsiEta[0]) * (1 - wspolrzednePC4KsiEta[1]);
        PC4NToTabC[0][2] = 0.25 * (1 + wspolrzednePC4KsiEta[0]) * (1 + wspolrzednePC4KsiEta[1]);
        PC4NToTabC[0][3] = 0.25 * (1 - wspolrzednePC4KsiEta[0]) * (1 + wspolrzednePC4KsiEta[1]);
    }
};



// definicje funkcji do calkowania
double fun1(double x) {
    return 5 * x * x + 3 * x + 6;
}

double fun2(double x, double y) {
    return 5 * x * x * y * y + 3 * x * y + 6;
}


// dla schematu 2punktowego
double weight2p[2] = { 1, 1 }; // tablica wag dla 2punktowego
double tabNodes2p[2] = { -1 / sqrt(3), 1 / sqrt(3) }; // tablica wezlow dla 2punktowego

// dla schematu 3punktowego
double weight3p[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 }; // {w1, w2, w3} tablica wag dla 3punktowego
double tabNodes3p[3] = { -sqrt(3.0 / 5.0), 0, sqrt(3.0 / 5.0) }; // tablica wezlow dla 3punktowego


double Gauss1D(int ile_punktowy) {
    double suma = 0;
    switch (ile_punktowy) {
    case 2:
        for (int i = 0; i < ile_punktowy; i++)
            suma += weight2p[i] * fun1(tabNodes2p[i]); // zliczane są iloczyny wagi i wartości funkcji w tym punkcie 
        break;
    case 3:
        for (int i = 0; i < ile_punktowy; i++)
            suma += weight3p[i] * fun1(tabNodes3p[i]);
        break;
    default:
        break;
    }
    return suma;
}


// wynikiem całkowania w każdym z punktów całkowania jest prostopadłościan o polu podstawy definiowanym wartością wag. 
// waga razy waga definiuje pole podstawy a wysokość to wartość funkcji w punkcie całkowania.
double Gauss2D(int ile_punktowy) {
    double suma = 0;
    switch (ile_punktowy)
    {
    case 2:
        for (int i = 0; i < ile_punktowy; i++)
        {
            for (int j = 0; j < ile_punktowy; j++)
            {
                suma += weight2p[j] * weight2p[i] * fun2(tabNodes2p[i], tabNodes2p[j]); // waga * waga * wartość funkcji w punkcie
            }
        }
        break;
    case 3:
        for (int i = 0; i < ile_punktowy; i++)
        {
            for (int j = 0; j < ile_punktowy; j++)
            {
                suma += weight3p[j] * weight3p[i] * fun2(tabNodes3p[i], tabNodes3p[j]);
            }
        }
        break;
    default:
        break;
    }
    return suma;
}


// Obliczamy pochodne funkcji kształtu względem ksi oraz eta:
// Za pomocą funkcji kształtu wyznacza się wsp. pkt. całk. w układzie globalnym
// Poniżej wzory dla elementu 4 węzłowego
double pochodnaKsi(int i, point p) {
    if (i == 0)
        return -1. / 4 * (1. - p.eta);
    if (i == 1)
        return 1. / 4 * (1. - p.eta);
    if (i == 2)
        return 1. / 4 * (1. + p.eta);
    if (i == 3)
        return -1. / 4 * (1. + p.eta);
}

double pochodnaEta(int i, point p) {
    if (i == 0)
        return (-1. / 4) * (1. - p.ksi);
    if (i == 1)
        return (-1. / 4) * (1. + p.ksi);
    if (i == 2)
        return (1. / 4) * (1. + p.ksi);
    if (i == 3)
        return (1. / 4) * (1. - p.ksi);
}


void showDetTab(double tab[])
{
    for (int i = 0; i < 4; i++)
    {
        if (i == 2)
        {
            cout << endl;
        }
        cout << "[" << tab[i] << "] ";
    }
    cout << endl << endl;
}


void showTab4(double tab[])
{
    for (int i = 0; i < 4; i++)
    {
        cout << "[" << tab[i] << "] " << endl;
    }
}

void showTab4x4(double tab[4][4])
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << "[" << tab[i][j] << "] ";
        }
        cout << endl;
    }
    cout << endl;
}

void showTab4x4Dynamic(double** tab)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << "[" << tab[i][j] << "] ";
        }
        cout << endl;
    }
    cout << endl;
}

void showTab2x2(double tab[2][2])
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            cout << "[" << tab[i][j] << "] ";
        }
        cout << endl;
    }
    cout << endl;
}

void showTab2x2Dynamic(double** tab)
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            cout << "[" << tab[i][j] << "] ";
        }
        cout << endl;
    }
    cout << endl;
}

void showTab2x4(double tab[2][4])
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << "[" << tab[i][j] << "] ";
        }
        cout << endl;
    }
}

void showTab1x4(double tab[1][4])
{
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << "[" << tab[i][j] << "] ";
        }
        cout << endl;
    }
}

void showTab4Dynamic(double* tab)
{
    for (int i = 0; i < 4; i++)
    {
        cout << "[" << tab[i] << "] " << endl;
    }
    cout << endl;
}

void showTab4x2(double tab[4][2])
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            cout << "[" << tab[i][j] << "] ";
        }
        cout << endl;
    }
}

void showTab16x16(double tab[16][16])
{
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            cout << "[" << tab[i][j] << "] ";
        }
        cout << endl;
    }
    cout << endl;
}

void showTab16x16Dynamic(double** tab)
{
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            cout << "[" << tab[i][j] << "] ";
        }
        cout << endl;
    }
    cout << endl;

}



// transoponowanie tablicy 1x4
double  transponujTablice1x4(double tab[1][4], double tabTransponowana[4][1])
{
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            tabTransponowana[j][i] = tab[i][j];
        }
    }
    return tabTransponowana[4][1];
}

// transoponowanie tablicy 4x4
double  transponujTablice4x4(double tab[4][4], double tabTransponowana[4][4])
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            tabTransponowana[j][i] = tab[i][j];
        }
    }
    return tabTransponowana[4][4];
}

// transponowanie tablicy 2x4
double transponujTablice2x4(double tab[2][4], double tabTransponowana[4][2])
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            tabTransponowana[j][i] = tab[i][j];
        }
    }
    return tabTransponowana[4][2];
}


void Jacobian(UniversalElement universalElement, element4 element, grid siatka, GlobalData globalData, double* macierzJakobianuPC1, double* macierzJakobianuPC2, double* macierzJakobianuPC3, double* macierzJakobianuPC4)
{
    double dXdKsiTab[4];
    double dXdEtaTab[4];
    double dYdKsiTab[4];
    double dYdEtaTab[4];

    double odwrotnoscDetJTab[4];

    int j = 0;
    for (int i = 0; i < 4; i++)
    {
        dXdKsiTab[i] = element.ksi[i][0] * siatka.tabNodes[siatka.tabElements[i].ID[0] - 1].y + element.ksi[i][1] * siatka.tabNodes[siatka.tabElements[i].ID[1] - 1].y + element.ksi[i][2] * siatka.tabNodes[siatka.tabElements[i].ID[2] - 1].y + element.ksi[i][3] * siatka.tabNodes[siatka.tabElements[i].ID[3] - 1].y;
        dXdEtaTab[i] = element.eta[i][0] * siatka.tabNodes[siatka.tabElements[i].ID[0] - 1].x + element.eta[i][1] * siatka.tabNodes[siatka.tabElements[i].ID[1] - 1].x + element.eta[i][2] * siatka.tabNodes[siatka.tabElements[i].ID[2] - 1].x + element.eta[i][3] * siatka.tabNodes[siatka.tabElements[i].ID[3] - 1].x;
        dYdKsiTab[i] = element.ksi[i][0] * siatka.tabNodes[siatka.tabElements[i].ID[0] - 1].x + element.ksi[i][1] * siatka.tabNodes[siatka.tabElements[i].ID[1] - 1].x + element.ksi[i][2] * siatka.tabNodes[siatka.tabElements[i].ID[2] - 1].x + element.ksi[i][3] * siatka.tabNodes[siatka.tabElements[i].ID[3] - 1].x;
        dYdEtaTab[i] = element.eta[i][0] * siatka.tabNodes[siatka.tabElements[i].ID[0] - 1].y + element.eta[i][1] * siatka.tabNodes[siatka.tabElements[i].ID[1] - 1].y + element.eta[i][2] * siatka.tabNodes[siatka.tabElements[i].ID[2] - 1].y + element.eta[i][3] * siatka.tabNodes[siatka.tabElements[i].ID[3] - 1].y;

        // Różnica długości wag związana jest z Jakobianem przekształcenia 1D (detJ)
        detJTab[i] = dYdEtaTab[i] * dYdKsiTab[i] - dXdKsiTab[i] * dXdEtaTab[i];
        //cout << endl << "det przy jakobianie:" << detJTab[i] << endl;
        //cout << endl << "dXdKsiTab:" << dXdKsiTab[i] << endl;
        //cout << endl << "dXdEtaTab:" << dXdEtaTab[i] << endl;
        //cout << endl << "dYdKsiTab:" << dYdKsiTab[i] << endl;
        //cout << endl << "dYdEtaTab:" << dYdEtaTab[i] << endl;

        odwrotnoscDetJTab[i] = 1 / detJTab[i];

        //cout << "jakobiany" << endl; 
        macierzJakobianuPC1[i] = odwrotnoscDetJTab[i] * dYdEtaTab[i];
        macierzJakobianuPC2[i] = odwrotnoscDetJTab[i] * dXdKsiTab[i];
        macierzJakobianuPC3[i] = odwrotnoscDetJTab[i] * dXdEtaTab[i];
        macierzJakobianuPC4[i] = odwrotnoscDetJTab[i] * dYdKsiTab[i];
        //cout << macierzJakobianuPC1[i] << " " << macierzJakobianuPC2[i] << " \n" << macierzJakobianuPC3[i] << " " << macierzJakobianuPC4[i] << endl;
    }
}


double** calculateTabH(pcKsiEta pcKsiEta, grid siatka, GlobalData globalData, UniversalElement universalElement, element4 element, double* macierzJakobianuPC1, double* macierzJakobianuPC2, double* macierzJakobianuPC3, double* macierzJakobianuPC4)
{
    double HPC1[4][4];
    double HPC2[4][4];
    double HPC3[4][4];
    double HPC4[4][4];
    double** H = new double* [siatka.nN];
    for (int i = 0; i < siatka.nN; i++)
    {
        H[i] = new double[siatka.nN];
    }

    double iloczynXTabPC1[4][4];
    double iloczynYTabPC1[4][4];
    double iloczynXTabPC2[4][4];
    double iloczynYTabPC2[4][4];
    double iloczynXTabPC3[4][4];
    double iloczynYTabPC3[4][4];
    double iloczynXTabPC4[4][4];
    double iloczynYTabPC4[4][4];

    double dNdXTab[4][4];
    double dNdYTab[4][4];

    double dNdXTabTransponowana[4][4];
    double dNdYTabTransponowana[4][4];

    //cout << "macierzJakobianuPC1" << endl;
    //showTab4(macierzJakobianuPC1);
    //cout << "macierzJakobianuPC2" << endl;
    //showTab4(macierzJakobianuPC2);
    //cout << "macierzJakobianuPC3" << endl;
    //showTab4(macierzJakobianuPC3);
    //cout << "macierzJakobianuPC4" << endl;
    //showTab4(macierzJakobianuPC4);


    // Tablica dNdX
    //cout << endl << "Tablica dNdX:" << endl;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            double dNdX = macierzJakobianuPC1[i] * element.ksi[i][j] + macierzJakobianuPC2[i] * element.eta[i][j];
            dNdXTab[i][j] = dNdX;
        }
    }

    // Tablica dNdY
   // cout << endl << "Tablica dNdY:" << endl;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            double dNdY = macierzJakobianuPC3[i] * element.ksi[i][j] + macierzJakobianuPC4[i] * element.eta[i][j];
            dNdYTab[i][j] = dNdY;
        }
    }

    transponujTablice4x4(dNdXTab, dNdXTabTransponowana);
    transponujTablice4x4(dNdYTab, dNdYTabTransponowana);



    for (int i = 0; i < globalData.nPC; i++) {
        for (int j = 0; j < globalData.nPC; j++) {
            iloczynXTabPC1[i][j] = dNdXTab[0][i] * dNdXTabTransponowana[j][0];
            iloczynYTabPC1[i][j] = dNdYTab[0][i] * dNdYTabTransponowana[j][0];
            HPC1[i][j] = globalData.conductivity * (iloczynXTabPC1[i][j] + iloczynYTabPC1[i][j]) * detJTab[0];


            iloczynXTabPC2[i][j] = dNdXTab[1][i] * dNdXTabTransponowana[j][1];
            iloczynYTabPC2[i][j] = dNdYTab[1][i] * dNdYTabTransponowana[j][1];
            HPC2[i][j] = globalData.conductivity * (iloczynXTabPC2[i][j] + iloczynYTabPC2[i][j]) * detJTab[1];

            iloczynXTabPC3[i][j] = dNdXTab[2][i] * dNdXTabTransponowana[j][2];
            iloczynYTabPC3[i][j] = dNdYTab[2][i] * dNdYTabTransponowana[j][2];
            HPC3[i][j] = globalData.conductivity * (iloczynXTabPC3[i][j] + iloczynYTabPC3[i][j]) * detJTab[2];

            iloczynXTabPC4[i][j] = dNdXTab[3][i] * dNdXTabTransponowana[j][3];
            iloczynYTabPC4[i][j] = dNdYTab[3][i] * dNdYTabTransponowana[j][3];
            HPC4[i][j] = globalData.conductivity * (iloczynXTabPC4[i][j] + iloczynYTabPC4[i][j]) * detJTab[3];

            H[i][j] = HPC1[i][j] + HPC2[i][j] + HPC3[i][j] + HPC4[i][j];
        }
    }
    //cout << "Macierz H (zsumowane macierze H lokalne):" << endl;
    //showTab4x4Dynamic(H);

    return H;
}



double* wektorP(UniversalElement universalElement, pcKsiEta pcKsiEta, grid siatka, GlobalData globalData, double*** HBC, double** PDwuWym)
{
    double wspolrzednePC1KsiEta[2][2] = { {-pcKsiEta.b, -pcKsiEta.a}, {pcKsiEta.b, -pcKsiEta.a} };
    double wspolrzednePC2KsiEta[2][2] = { {pcKsiEta.a, pcKsiEta.b}, {pcKsiEta.a, -pcKsiEta.b} };
    double wspolrzednePC3KsiEta[2][2] = { {-pcKsiEta.b, pcKsiEta.a}, {pcKsiEta.b, pcKsiEta.a} };
    double wspolrzednePC4KsiEta[2][2] = { {-pcKsiEta.a, pcKsiEta.b}, {-pcKsiEta.a, -pcKsiEta.b} };

    double PC1N[2][4];
    double PC2N[2][4];
    double PC3N[2][4];
    double PC4N[2][4];

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            double N1_PC4 = 0.25 * (1 - wspolrzednePC4KsiEta[i][0]) * (1 - wspolrzednePC4KsiEta[i][1]);
            double N2_PC4 = 0.25 * (1 + wspolrzednePC4KsiEta[i][0]) * (1 - wspolrzednePC4KsiEta[i][1]);
            double N3_PC4 = 0.25 * (1 + wspolrzednePC4KsiEta[i][0]) * (1 + wspolrzednePC4KsiEta[i][1]);
            double N4_PC4 = 0.25 * (1 - wspolrzednePC4KsiEta[i][0]) * (1 + wspolrzednePC4KsiEta[i][1]);

            double N1_PC1 = 0.25 * (1 - wspolrzednePC1KsiEta[i][0]) * (1 - wspolrzednePC1KsiEta[i][1]);
            double N2_PC1 = 0.25 * (1 + wspolrzednePC1KsiEta[i][0]) * (1 - wspolrzednePC1KsiEta[i][1]);
            double N3_PC1 = 0.25 * (1 + wspolrzednePC1KsiEta[i][0]) * (1 + wspolrzednePC1KsiEta[i][1]);
            double N4_PC1 = 0.25 * (1 - wspolrzednePC1KsiEta[i][0]) * (1 + wspolrzednePC1KsiEta[i][1]);

            double N1_PC2 = 0.25 * (1 - wspolrzednePC2KsiEta[i][0]) * (1 - wspolrzednePC2KsiEta[i][1]);
            double N2_PC2 = 0.25 * (1 + wspolrzednePC2KsiEta[i][0]) * (1 - wspolrzednePC2KsiEta[i][1]);
            double N3_PC2 = 0.25 * (1 + wspolrzednePC2KsiEta[i][0]) * (1 + wspolrzednePC2KsiEta[i][1]);
            double N4_PC2 = 0.25 * (1 - wspolrzednePC2KsiEta[i][0]) * (1 + wspolrzednePC2KsiEta[i][1]);

            double N1_PC3 = 0.25 * (1 - wspolrzednePC3KsiEta[i][0]) * (1 - wspolrzednePC3KsiEta[i][1]);
            double N2_PC3 = 0.25 * (1 + wspolrzednePC3KsiEta[i][0]) * (1 - wspolrzednePC3KsiEta[i][1]);
            double N3_PC3 = 0.25 * (1 + wspolrzednePC3KsiEta[i][0]) * (1 + wspolrzednePC3KsiEta[i][1]);
            double N4_PC3 = 0.25 * (1 - wspolrzednePC3KsiEta[i][0]) * (1 + wspolrzednePC3KsiEta[i][1]);
            PC4N[i][0] = N1_PC4;
            PC4N[i][1] = N2_PC4;
            PC4N[i][2] = N3_PC4;
            PC4N[i][3] = N4_PC4;

            PC1N[i][0] = N1_PC1;
            PC1N[i][1] = N2_PC1;
            PC1N[i][2] = N3_PC1;
            PC1N[i][3] = N4_PC1;

            PC2N[i][0] = N1_PC2;
            PC2N[i][1] = N2_PC2;
            PC2N[i][2] = N3_PC2;
            PC2N[i][3] = N4_PC2;

            PC3N[i][0] = N1_PC3;
            PC3N[i][1] = N2_PC3;
            PC3N[i][2] = N3_PC3;
            PC3N[i][3] = N4_PC3;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //obliczanie HBC

    double PC1NTransponowana[4][2];
    double PC2NTransponowana[4][2];
    double PC3NTransponowana[4][2];
    double PC4NTransponowana[4][2];

    // transponowanie tablicy N
    // dla PC1
    transponujTablice2x4(PC1N, PC1NTransponowana);
    // dla PC2
    transponujTablice2x4(PC2N, PC2NTransponowana);
    //dla PC3
    transponujTablice2x4(PC3N, PC3NTransponowana);
    // dla PC4
    transponujTablice2x4(PC4N, PC4NTransponowana);

    double iloczynPC1N_1[4][4];
    double iloczynPC1N_2[4][4];
    double iloczynPC2N_1[4][4];
    double iloczynPC2N_2[4][4];
    double iloczynPC3N_1[4][4];
    double iloczynPC3N_2[4][4];
    double iloczynPC4N_1[4][4];
    double iloczynPC4N_2[4][4];

    double PC4iloczynPart1[4][4];
    double PC4iloczynPart2[4][4];
    double PC1iloczynPart1[4][4];
    double PC1iloczynPart2[4][4];
    double PC2iloczynPart1[4][4];
    double PC2iloczynPart2[4][4];
    double PC3iloczynPart1[4][4];
    double PC3iloczynPart2[4][4];

    double PC4PiloczynPart1[4][4];
    double PC4PiloczynPart2[4][4];
    double PC1PiloczynPart1[4][4];
    double PC1PiloczynPart2[4][4];
    double PC2PiloczynPart1[4][4];
    double PC2PiloczynPart2[4][4];
    double PC3PiloczynPart1[4][4];
    double PC3PiloczynPart2[4][4];

    double HBC1[4][4];
    double HBC2[4][4];
    double HBC3[4][4];
    double HBC4[4][4];
    double HBCElement[4][4];
    double HBCAll[4][4];

    double detJ_HBC = siatka.tabNodes[siatka.tabElements[0].ID[3] - 1].y / 2; // jakobian dla HBC


    double P1[4];
    double P2[4];
    double P3[4];
    double P4[4];
    double PElement[4];
    double* P = new double[siatka.nN];


    double emptyArray[4][4];
    double emptyArrayP[4];
    for (int i = 0; i < 4; i++) {
        emptyArrayP[i] = 0;
        PElement[i] = 0;
        for (int j = 0; j < 4; j++) {
            emptyArray[i][j] = 0;
            HBC1[i][j] = 0;
            HBC2[i][j] = 0;
            HBC3[i][j] = 0;
            HBC4[i][j] = 0;
            HBCElement[i][j] = 0;
            HBCAll[i][j] = 0;
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {


            // HBC dla PC1
            iloczynPC4N_1[i][j] = PC4NTransponowana[i][0] * PC4N[0][j];
            PC4iloczynPart1[i][j] = weight2p[0] * globalData.alfa * iloczynPC4N_1[i][j];

            iloczynPC4N_2[i][j] = PC4NTransponowana[i][1] * PC4N[1][j];
            PC4iloczynPart2[i][j] = weight2p[0] * globalData.alfa * iloczynPC4N_2[i][j];
            HBC1[i][j] = (PC4iloczynPart1[i][j] + PC4iloczynPart2[i][j]) * detJ_HBC;


            //wektor P dla PC1
            PC4PiloczynPart1[i][j] = weight2p[0] * globalData.alfa * PC4NTransponowana[i][0] * globalData.t0;
            PC4PiloczynPart2[i][j] = weight2p[0] * globalData.alfa * PC4NTransponowana[i][1] * globalData.t0;
            P1[i] = (PC4PiloczynPart1[i][j] + PC4PiloczynPart2[i][j]) * detJ_HBC;


            // HBC dla PC2
            iloczynPC1N_1[i][j] = PC1NTransponowana[i][0] * PC1N[0][j];
            PC1iloczynPart1[i][j] = weight2p[0] * globalData.alfa * iloczynPC1N_1[i][j];

            iloczynPC1N_2[i][j] = PC1NTransponowana[i][1] * PC1N[1][j];
            PC1iloczynPart2[i][j] = weight2p[0] * globalData.alfa * iloczynPC1N_2[i][j];
            HBC2[i][j] = (PC1iloczynPart1[i][j] + PC1iloczynPart2[i][j]) * detJ_HBC;

            //wektor P dla PC2
            PC1PiloczynPart1[i][j] = weight2p[0] * globalData.alfa * PC1NTransponowana[i][0] * globalData.t0;
            PC1PiloczynPart2[i][j] = weight2p[0] * globalData.alfa * PC1NTransponowana[i][1] * globalData.t0;
            P2[i] = (PC1PiloczynPart1[i][j] + PC1PiloczynPart2[i][j]) * detJ_HBC;

            // HBC dla PC3
            iloczynPC2N_1[i][j] = PC2NTransponowana[i][0] * PC2N[0][j];
            PC2iloczynPart1[i][j] = weight2p[0] * globalData.alfa * iloczynPC2N_1[i][j];

            iloczynPC2N_2[i][j] = PC2NTransponowana[i][1] * PC2N[1][j];
            PC2iloczynPart2[i][j] = weight2p[0] * globalData.alfa * iloczynPC2N_2[i][j];
            HBC3[i][j] = (PC2iloczynPart1[i][j] + PC2iloczynPart2[i][j]) * detJ_HBC;

            //wektor P dla PC3
            PC2PiloczynPart1[i][j] = weight2p[0] * globalData.alfa * PC2NTransponowana[i][0] * globalData.t0;
            PC2PiloczynPart2[i][j] = weight2p[0] * globalData.alfa * PC2NTransponowana[i][1] * globalData.t0;
            P3[i] = (PC2PiloczynPart1[i][j] + PC2PiloczynPart2[i][j]) * detJ_HBC;

            // HBC dla PC4
            iloczynPC3N_1[i][j] = PC3NTransponowana[i][0] * PC3N[0][j];
            PC3iloczynPart1[i][j] = weight2p[0] * globalData.alfa * iloczynPC3N_1[i][j];

            iloczynPC3N_2[i][j] = PC3NTransponowana[i][1] * PC3N[1][j];
            PC3iloczynPart2[i][j] = weight2p[0] * globalData.alfa * iloczynPC3N_2[i][j];
            HBC4[i][j] = (PC3iloczynPart1[i][j] + PC3iloczynPart2[i][j]) * detJ_HBC;


            //wektor P dla PC4
            PC3PiloczynPart1[i][j] = weight2p[0] * globalData.alfa * PC3NTransponowana[i][0] * globalData.t0;
            PC3PiloczynPart2[i][j] = weight2p[0] * globalData.alfa * PC3NTransponowana[i][1] * globalData.t0;
            P4[i] = (PC3PiloczynPart1[i][j] + PC3PiloczynPart2[i][j]) * detJ_HBC;

            P[i] = P1[i] + P2[i] + P3[i] + P4[i];
        }
    }



    for (int e = 0; e < siatka.nE; e++)
    {


        if (siatka.tabNodes[siatka.tabElements[e].ID[0] - 1].y == 0 && siatka.tabNodes[siatka.tabElements[e].ID[0] - 1].x == 0) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    HBCElement[i][j] = HBC1[i][j] + HBC2[i][j];
                    PElement[i] = P1[i] + P2[i];
                }
            }
        }
        else if (siatka.tabNodes[siatka.tabElements[e].ID[1] - 1].y == 0 && siatka.tabNodes[siatka.tabElements[e].ID[1] - 1].x == siatka.B) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    HBCElement[i][j] = HBC2[i][j] + HBC3[i][j];
                    PElement[i] = P2[i] + P3[i];
                }
            }
        }
        else if (siatka.tabNodes[siatka.tabElements[e].ID[2] - 1].y == siatka.H && siatka.tabNodes[siatka.tabElements[e].ID[2] - 1].x == siatka.B) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    HBCElement[i][j] = HBC3[i][j] + HBC4[i][j];
                    PElement[i] = P3[i] + P4[i];

                }
            }
        }
        else if (siatka.tabNodes[siatka.tabElements[e].ID[3] - 1].y == siatka.H && siatka.tabNodes[siatka.tabElements[e].ID[3] - 1].x == 0) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    HBCElement[i][j] = HBC1[i][j] + HBC4[i][j];
                    PElement[i] = P1[i] + P4[i];
                }
            }
        }

        else if (siatka.tabNodes[siatka.tabElements[e].ID[0] - 1].x == 0 &&
            siatka.tabNodes[siatka.tabElements[e].ID[0] - 1].y > 0 &&
            siatka.tabNodes[siatka.tabElements[e].ID[3] - 1].y < siatka.H) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    HBCElement[i][j] = HBC1[i][j] + emptyArray[i][j];
                    PElement[i] = P1[i] + emptyArrayP[i];

                }
            }
        }
        else if (siatka.tabNodes[siatka.tabElements[e].ID[1] - 1].x == siatka.B &&
            siatka.tabNodes[siatka.tabElements[e].ID[1] - 1].y > 0 &&
            siatka.tabNodes[siatka.tabElements[e].ID[2] - 1].y < siatka.H) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    HBCElement[i][j] = HBC3[i][j] + emptyArray[i][j];
                    PElement[i] = P3[i] + emptyArrayP[i];

                }
            }
        }
        else if (siatka.tabNodes[siatka.tabElements[e].ID[0] - 1].y == 0 &&
            siatka.tabNodes[siatka.tabElements[e].ID[0] - 1].x > 0 &&
            siatka.tabNodes[siatka.tabElements[e].ID[1] - 1].x < siatka.B) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    HBCElement[i][j] = HBC2[i][j] + emptyArray[i][j];
                    PElement[i] = P2[i] + emptyArrayP[i];

                }
            }
        }
        else if (siatka.tabNodes[siatka.tabElements[e].ID[3] - 1].y == siatka.H &&
            siatka.tabNodes[siatka.tabElements[e].ID[3] - 1].x > 0 &&
            siatka.tabNodes[siatka.tabElements[e].ID[2] - 1].x < siatka.B) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    HBCElement[i][j] = HBC4[i][j] + emptyArray[i][j];
                    PElement[i] = P4[i] + emptyArrayP[i];
                }
            }
        }

        // showTab4x4(HBCElement);
         //showTab4(PElement);



        for (int i = 0; i < 4; i++) {
            PDwuWym[e][i] = PElement[i];
            for (int j = 0; j < 4; j++) {
                HBC[e][i][j] = HBCElement[i][j];
            }
        }

        for (int i = 0; i < 4; i++) {
            PElement[i] = 0;
            for (int j = 0; j < 4; j++) {
                HBCElement[i][j] = 0;
            }
        }

    }



    // GENERYCZNOSC MACIERZY P //

    // tworzenie wektora globalnego P, rozmiar tablicy = liczbawezlow 
    double* GlobalnyWektorP = new double[siatka.nN];
    for (int i = 0; i < siatka.nN; i++)
    {
        GlobalnyWektorP[i] = 0;
    }


    return GlobalnyWektorP;
}


// generyczność wektora P
double* generycznoscP(grid siatka, double* Pglobal, double** P) // argumenty: siatka, pusta tablica na wyniki , tablica P zawierająca tablice P poszczególnych elementów
{
    for (int e = 0; e < siatka.nE; e++)
    {
        for (int i = 0; i < 4; i++) {
            Pglobal[siatka.tabElements[e].ID[i] - 1] += P[e][i];
        }
    }

    return Pglobal;
}


double** calculateTabC(pcKsiEta pcKsiEta, grid siatka, GlobalData globalData, UniversalElement universalElement)
{
    double PC1NToTabCTransponowana[4][1];
    double PC2NToTabCTransponowana[4][1];
    double PC3NToTabCTransponowana[4][1];
    double PC4NToTabCTransponowana[4][1];
    double C1[4][4];
    double C2[4][4];
    double C3[4][4];
    double C4[4][4];

    double** C = new double* [siatka.nN];
    for (int i = 0; i < siatka.nN; i++)
    {
        C[i] = new double[siatka.nN];
    }

    transponujTablice1x4(universalElement.PC1NToTabC, PC1NToTabCTransponowana);
    transponujTablice1x4(universalElement.PC2NToTabC, PC2NToTabCTransponowana);
    transponujTablice1x4(universalElement.PC3NToTabC, PC3NToTabCTransponowana);
    transponujTablice1x4(universalElement.PC4NToTabC, PC4NToTabCTransponowana);


    for (int i = 0; i < globalData.nPC; i++)
    {
        for (int j = 0; j < globalData.nPC; j++)
        {
            C1[i][j] = globalData.c * globalData.ro * (PC1NToTabCTransponowana[i][0] * universalElement.PC1NToTabC[0][j]) * detJTab[0];
            C2[i][j] = globalData.c * globalData.ro * (PC2NToTabCTransponowana[i][0] * universalElement.PC2NToTabC[0][j]) * detJTab[0];
            C3[i][j] = globalData.c * globalData.ro * (PC3NToTabCTransponowana[i][0] * universalElement.PC3NToTabC[0][j]) * detJTab[0];
            C4[i][j] = globalData.c * globalData.ro * (PC4NToTabCTransponowana[i][0] * universalElement.PC4NToTabC[0][j]) * detJTab[0];
            C[i][j] = C1[i][j] + C2[i][j] + C3[i][j] + C4[i][j];
        }
    }

    //cout << endl << "Macierz C dla PC1:" << endl;
    //showTab4x4(C1);
    //cout << endl << "Macierz C dla PC2:" << endl;
    //showTab4x4(C2);
    //cout << endl << "Macierz C dla PC3:" << endl;
    //showTab4x4(C3);
    //cout << endl << "Macierz C dla PC4:" << endl;
    //showTab4x4(C4);
    //cout << endl << "Macierz C:" << endl;
    //showTab4x4Dynamic(C);

    return C;
}

// funkcja agregująca tablice przesłaną jako argument
double** agregacja(double** tablicaDoAgregacji, grid siatka)
{
    // tworzenie tablicy globalnej, rozmiar tablicy = liczbawezlow x liczbawezlow
    double** TabGlobalnaAgregacja = new double* [siatka.nN];
    for (int i = 0; i < siatka.nN; i++)
    {
        TabGlobalnaAgregacja[i] = new double[siatka.nN];
    }

    for (int i = 0; i < siatka.nN; i++) {
        for (int j = 0; j < siatka.nN; j++) {
            TabGlobalnaAgregacja[i][j] = 0;
        }
    }

    for (int e = 0; e < siatka.nE; e++)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                // Tab[1-1][1-1]
                TabGlobalnaAgregacja[siatka.tabElements[e].ID[i] - 1][siatka.tabElements[e].ID[j] - 1] += tablicaDoAgregacji[i][j];
            }
        }
    }

    return TabGlobalnaAgregacja;
}


// agregacja do HBC, jako argument przyjmuje tablice tablic HBC dla poszczególnych elementów
double** agregacjaHBC(double*** tablicaDoAgregacji, grid siatka)
{
    // tworzenie tablicy globalnej, rozmiar tablicy = liczbawezlow x liczbawezlow
    double** TabGlobalnaAgregacja = new double* [siatka.nN];
    for (int i = 0; i < siatka.nN; i++)
    {
        TabGlobalnaAgregacja[i] = new double[siatka.nN];
    }

    for (int i = 0; i < siatka.nN; i++) {
        for (int j = 0; j < siatka.nN; j++) {
            TabGlobalnaAgregacja[i][j] = 0;
        }
    }

    for (int e = 0; e < siatka.nE; e++)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                // Tab[1-1][1-1]
                TabGlobalnaAgregacja[siatka.tabElements[e].ID[i] - 1][siatka.tabElements[e].ID[j] - 1] += tablicaDoAgregacji[e][i][j];
            }
        }
    }

    return TabGlobalnaAgregacja;
}


void calculateMatrixH(GlobalData globalData, grid& siatka, double** C, double** H, double** matrixH, double** tabGlobalH, double** HBC)
{
    for (int i = 0; i < siatka.nN; i++)
    {
        for (int j = 0; j < siatka.nN; j++)
        {
            tabGlobalH[i][j] = HBC[i][j] + H[i][j];
        }
    }


    for (int i = 0; i < siatka.nN; i++)
    {
        for (int j = 0; j < siatka.nN; j++)
        {
            matrixH[i][j] += tabGlobalH[i][j] + (C[i][j] / globalData.sst); // Hglobal = HBC + H
        }
    }

    cout << endl << "Ostatateczny Matrix H:" << endl;
    for (int i = 0; i < siatka.nN; i++) {
        for (int j = 0; j < siatka.nN; j++) {
            std::cout << std::setprecision(4) << "[" << matrixH[i][j] << "]";
        }cout << endl;
    }
}


void P_vector(GlobalData globalData, grid& siatka, double** CdoRownania, double* PDoRownania, double* vectorP)
{
    for (int i = 0; i < siatka.nN; i++)
    {
        vectorP[i] = 0;
    }
    for (int i = 0; i < siatka.nN; i++)
    {
        vectorP[i] += PDoRownania[i];

        for (int j = 0; j < siatka.nN; j++)
        {
            vectorP[i] += (CdoRownania[i][j] * siatka.tabNodes[j].temp0) / globalData.sst;
        }
    }

    cout << endl << "P_vector:" << endl;
    for (int i = 0; i < siatka.nN; i++)
    {
        cout << "[" << vectorP[i] << "] " << endl;
    }
    cout << endl;

}


void calculateHP(grid siatka, double** HPTab, double** matrixH, double* vectorP) {
    for (int i = 0; i < siatka.nN; i++)
    {

        for (int j = 0; j < siatka.nN + 1; j++)
        {
            if (j == siatka.nN)
            {
                HPTab[i][j] = vectorP[i];
            }
            else
                HPTab[i][j] = matrixH[i][j];
        }
    }

    cout << endl << "Macierz HP:" << endl;
    for (int i = 0; i < siatka.nN; i++)
    {
        for (int j = 0; j < siatka.nN; j++)
        {
            std::cout << std::setprecision(4) << "[" << HPTab[i][j] << "]";
        }
        cout << endl;
    }

}



//obliczanie ukladow rownan (HP), wynik w postaci wektora temperatur
double* EliminacjaGaussa(double** Tab, int size)
{
    double mnoznik = 0.0f;
    for (int k = 0; k < (size - 1); k++)
    {
        for (int i = (k + 1); i < size; i++)
        {
            if (Tab[k][k] == 0) // sprawdzenie czy na przekątnej nie ma 0
            {
                cout << "Na przekatnej macierzy wystapila wartosc 0 !!!" << endl;
                //return;
            }
            mnoznik = Tab[i][k] / Tab[k][k]; // obliczamy mnoznik
            for (int j = k; j < (size + 1); j++)
            {
                Tab[i][j] -= (mnoznik * Tab[k][j]); // od i-tego wiersza odejmujemy wiersz zerowy pomnożony przez mnożnik 
            }
        }
    }

    double* ResultTab = new double[size];
    for (int i = 0; i < size; i++)
    {
        ResultTab[i] = 1;
    }

    double l = 0;
    for (int i = (size - 1); i >= 0; i--)
    {
        l = Tab[i][size];
        for (int j = (size - 1); j > i; j--)
        {
            l -= Tab[i][j] * ResultTab[j];
        }
        ResultTab[i] = l / Tab[i][i];
    }

    // Wypisanie wyników
    cout << endl << "Rozwiazanie ukladu rownan:" << endl;
    for (int i = 0; i < size; i++)
    {
        cout << "temp_" << i << " = " << ResultTab[i] << endl;
    }


    return ResultTab;




}



int main()
{
    GlobalData globalData; // instancja typu GlobalData ze zmiennymi
    UniversalElement universalElement(2.0); // utworzenie elementu uniwersalnego
    element e;
    //int nPC = 4; // liczba punktów całkowania 
    element4 element4;
    pcKsiEta pcKsiEta;


    double** HBC = new double* [4];

    for (int i = 0; i < 4; i++) {
        HBC[i] = new double[4];
    }

    // siatka mes
    cout << endl << "Siatka MES:" << endl;
    grid siatka;
    for (int i = 0; i < siatka.nE; i++)
    {
        cout << "Element nr." << i + 1 << ": " <<
            siatka.tabElements[i].ID[0] << " " <<
            siatka.tabElements[i].ID[1] << " " <<
            siatka.tabElements[i].ID[2] << " " <<
            siatka.tabElements[i].ID[3];
        cout << endl;
    }
    for (int i = 0; i < siatka.nN; i++)
    {
        cout << "Wezel nr." << i + 1 << ": " << "x: " << siatka.tabNodes[i].x << ", y: " << siatka.tabNodes[i].y << endl;
    }

    ////////////////////////// wyznaczanie warunkow brzegowych na scianach elementow /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < siatka.nE; i++)
    {
        // ELEMENTY NA ROGACH SIATKI - warunek brzegowy na dwoch scianach

        // warunek określający element w rogu na samym dole po lewej stronie, ponieważ id[0]: x=0, y=0
        if (siatka.tabNodes[siatka.tabElements[i].ID[0] - 1].y == 0 && siatka.tabNodes[siatka.tabElements[i].ID[0] - 1].x == 0)
        {
            // na dwoch scianach wystepuje warunek brzegowy poniewaz jest to element w rogu 
            siatka.tabElements[i].TabSciany[0].BC = true; // flage BC na scianie ustawiam na true ponieważ występuje na niej warunek brzegowy
            siatka.tabElements[i].TabSciany[3].BC = true;
        }
        // warunek określający element w rogu na samym dole po prawej stronie, ponieważ id[1]: x=g.B, y=0
        else if (siatka.tabNodes[siatka.tabElements[i].ID[1] - 1].y == 0 && siatka.tabNodes[siatka.tabElements[i].ID[1] - 1].x == siatka.B)
        {
            // na dwoch scianach wystepuje warunek brzegowy poniewaz jest to element w rogu 
            siatka.tabElements[i].TabSciany[0].BC = true; // flage BC na scianie ustawiam na true ponieważ występuje na niej warunek brzegowy
            siatka.tabElements[i].TabSciany[1].BC = true;
        }
        // warunek określający element w rogu na samej górze po prawej stronie, ponieważ id[2]: x=g.B, y=g.H
        else if (siatka.tabNodes[siatka.tabElements[i].ID[2] - 1].y == siatka.H && siatka.tabNodes[siatka.tabElements[i].ID[2] - 1].x == siatka.B)
        {
            // na dwoch scianach wystepuje warunek brzegowy poniewaz jest to element w rogu 
            siatka.tabElements[i].TabSciany[1].BC = true; // flage BC na scianie ustawiam na true ponieważ występuje na niej warunek brzegowy
            siatka.tabElements[i].TabSciany[2].BC = true;
        }
        // warunek określający element w rogu na samej górze po lewej stronie, ponieważ id[3]: x=0, y=g.H
        else if (siatka.tabNodes[siatka.tabElements[i].ID[3] - 1].y == siatka.H && siatka.tabNodes[siatka.tabElements[i].ID[3] - 1].x == 0)
        {
            // na dwoch scianach wystepuje warunek brzegowy poniewaz jest to element w rogu 
            siatka.tabElements[i].TabSciany[2].BC = true; // flage BC na scianie ustawiam na true ponieważ występuje na niej warunek brzegowy
            siatka.tabElements[i].TabSciany[3].BC = true;
        }

        // ELEMENTY NA SCIANACH SIATKI ALE NIE ROGI - - warunek brzegowy na jednej ścianie

        // warunek określający element na lewej scianie siatki bez naroznikow, ponieważ id[0]: x=0, y>0, id[3]: x<g.H
        else if (siatka.tabNodes[siatka.tabElements[i].ID[0] - 1].x == 0 && siatka.tabNodes[siatka.tabElements[i].ID[0] - 1].y > 0 && siatka.tabNodes[siatka.tabElements[i].ID[3] - 1].y < siatka.H)
        {
            // na jednej scianie wystepuje warunek brzegowy poniewaz jest to element na scianie ale nie w rogu 
            siatka.tabElements[i].TabSciany[3].BC = true; // flage BC na scianie ustawiam na true ponieważ występuje na niej warunek brzegowy
        }
        // warunek określający element na prawej scianie siatki bez naroznikow, ponieważ id[1]: x=g.B, y>0, id[2]: x<g.H
        else if (siatka.tabNodes[siatka.tabElements[i].ID[1] - 1].x == siatka.B && siatka.tabNodes[siatka.tabElements[i].ID[1] - 1].y > 0 && siatka.tabNodes[siatka.tabElements[i].ID[2] - 1].y < siatka.H)
        {
            // na jednej scianie wystepuje warunek brzegowy poniewaz jest to element na scianie ale nie w rogu 
            siatka.tabElements[i].TabSciany[1].BC = true; // flage BC na scianie ustawiam na true ponieważ występuje na niej warunek brzegowy
        }
        // warunek określający element na dolnej scianie siatki bez naroznikow, ponieważ id[0]: x>0, y=0, id[1]: x<g.B
        else if (siatka.tabNodes[siatka.tabElements[i].ID[0] - 1].y == 0 && siatka.tabNodes[siatka.tabElements[i].ID[0] - 1].x > 0 && siatka.tabNodes[siatka.tabElements[i].ID[1] - 1].x < siatka.B)
        {
            // na jednej scianie wystepuje warunek brzegowy poniewaz jest to element na scianie ale nie w rogu 
            siatka.tabElements[i].TabSciany[0].BC = true; // flage BC na scianie ustawiam na true ponieważ występuje na niej warunek brzegowy

        }
        // warunek określający element na górnej scianie siatki bez naroznikow, ponieważ id[3]: x>0, y=g.H, id[2]: x<g.B
        else if (siatka.tabNodes[siatka.tabElements[i].ID[3] - 1].y == siatka.H && siatka.tabNodes[siatka.tabElements[i].ID[3] - 1].x > 0 && siatka.tabNodes[siatka.tabElements[i].ID[2] - 1].x < siatka.B)
        {
            // na jednej scianie wystepuje warunek brzegowy poniewaz jest to element na scianie ale nie w rogu 
            siatka.tabElements[i].TabSciany[2].BC = true; // flage BC na scianie ustawiam na true ponieważ występuje na niej warunek brzegowy
        }
    }

    for (int i = 0; i < siatka.nE; i++)
    {
        cout << "BC dla element nr." << i + 1 << ": ";
        for (int j = 0; j < 4; j++)
        {
            if (siatka.tabElements[i].TabSciany[j].BC == true)
            {
                cout << "true, ";
            }
            else
            {
                cout << "false, ";
            }
        }
        cout << endl;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // całkowanie Gausa 
    cout << endl << "Calkowanie Gausa:" << endl;
    cout << "Gauss 1D dla dwoch punktow: " << Gauss1D(2) << endl;
    cout << "Gauss 1D dla trzech punktow: " << Gauss1D(3) << endl;
    cout << "Gauss 2D dla dwoch punktow: " << Gauss2D(2) << endl;
    cout << "Gauss 2D dla trzech punktow: " << Gauss2D(3) << endl;

    // Macierze pochodnych ksi, eta

    // punkty całkowania
    point pc1 = point(-1. / sqrt(3), -1. / sqrt(3));
    point pc2 = point(1. / sqrt(3), -1. / sqrt(3));
    point pc3 = point(1. / sqrt(3), 1. / sqrt(3));
    point pc4 = point(-1. / sqrt(3), 1. / sqrt(3));
    point p;

    //macierze pochodnych
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            switch (i) {
            case 0: p = pc1; break;
            case 1: p = pc2; break;
            case 2: p = pc3; break;
            case 3: p = pc4; break;
            default: break;
            }
            element4.ksi[i][j] = pochodnaKsi(j, p);
            element4.eta[i][j] = pochodnaEta(j, p);
        }
    //cout << endl << "Tablica pochodnych po ksi:" << endl;
    //showTab4x4Dynamic(element4.ksi);
    //cout << endl << "Tablica pochodnych po eta:" << endl;
    //showTab4x4Dynamic(element4.eta);

  /*  double* macierzJakobianuPC1 = new double[siatka.nN];
    double* macierzJakobianuPC2 = new double[siatka.nN];
    double* macierzJakobianuPC3 = new double[siatka.nN];
    double* macierzJakobianuPC4 = new double[siatka.nN];*/
    double* macierzJakobianuPC1 = new double[4];
    double* macierzJakobianuPC2 = new double[4];
    double* macierzJakobianuPC3 = new double[4];
    double* macierzJakobianuPC4 = new double[4];

    double** macierzJakobianuPC11 = new double* [2];
    double** macierzJakobianuPC22 = new double* [2];
    double** macierzJakobianuPC33 = new double* [2];
    double** macierzJakobianuPC44 = new double* [2];

    for (int i = 0; i < 2; i++)
    {
        macierzJakobianuPC11[i] = new double[2];
        macierzJakobianuPC22[i] = new double[2];
        macierzJakobianuPC33[i] = new double[2];
        macierzJakobianuPC44[i] = new double[2];

    }


    Jacobian(universalElement, element4, siatka, globalData, macierzJakobianuPC1, macierzJakobianuPC2, macierzJakobianuPC3, macierzJakobianuPC4);

    double* PDoRownania = new double[siatka.nN];
    double* vectorP = new double[siatka.nN];
    for (int i = 0; i < siatka.nN; i++)
    {
        vectorP[i] = 0;
    }

    //showTab1x4(universalElement.PC1NToTabC);
    //showTab1x4(universalElement.PC2NToTabC);
    //showTab1x4(universalElement.PC3NToTabC);
    //showTab1x4(universalElement.PC4NToTabC);

    double** C = new double* [siatka.nN];
    double** CdoRownania = new double* [siatka.nN];
    double** HBCdoRownania = new double* [siatka.nN];
    double** H = new double* [siatka.nN];
    double** HdoRownania = new double* [siatka.nN];
    double** tabGlobal1 = new double* [siatka.nN];
    double** tabGlobal3 = new double* [siatka.nN];
    double** matrixH = new double* [siatka.nN];
    double** tabGlobalH = new double* [siatka.nN];
    double** HPTab = new double* [siatka.nN];
    double** T2 = new double* [siatka.nN];
    for (int i = 0; i < siatka.nN; i++)
    {
        C[i] = new double[siatka.nN];
        CdoRownania[i] = new double[siatka.nN];
        HBCdoRownania[i] = new double[siatka.nN];
        HdoRownania[i] = new double[siatka.nN];
        H[i] = new double[siatka.nN];
        tabGlobal1[i] = new double[siatka.nN];
        tabGlobal3[i] = new double[siatka.nN];
        matrixH[i] = new double[siatka.nN];
        tabGlobalH[i] = new double[siatka.nN];
        HPTab[i] = new double[siatka.nN];
        T2[i] = new double[siatka.nN];
    }

    for (int i = 0; i < siatka.nN; i++)
    {
        for (int j = 0; j < siatka.nN; j++)
        {
            tabGlobal3[i][j] = 0;
            matrixH[i][j] = 0;
            tabGlobalH[i][j] = 0;
            HPTab[i][j] = 0;
            T2[i][j] = 0;
        }

    }

    C = calculateTabC(pcKsiEta, siatka, globalData, universalElement);
    CdoRownania = calculateTabC(pcKsiEta, siatka, globalData, universalElement);


    //cout << "Agregacja C:";
    //agregacja(C, g);
    CdoRownania = agregacja(C, siatka);


    H = calculateTabH(pcKsiEta, siatka, globalData, universalElement, element4, macierzJakobianuPC1, macierzJakobianuPC2, macierzJakobianuPC3, macierzJakobianuPC4);
    HdoRownania = agregacja(H, siatka);

    /*  cout << endl << "Agregacja macierzy H:" << endl;

      for (int i = 0; i < siatka.nN; i++) {
          for (int j = 0; j < siatka.nN; j++) {
              std::cout << std::setprecision(4) << "[" << HdoRownania[i][j] << "]";
          }cout << endl;
      }*/


    double*** HBCTrojWym = new double** [siatka.nE];


    for (int i = 0; i < siatka.nE; ++i) {
        HBCTrojWym[i] = new double* [4];
        for (int j = 0; j < 4; ++j) {
            HBCTrojWym[i][j] = new double[4];
        }
    }

    double** PDwuWym = new double* [siatka.nE];
    for (int i = 0; i < siatka.nE; i++)
    {
        PDwuWym[i] = new double[4];
    }

    PDoRownania = wektorP(universalElement, pcKsiEta, siatka, globalData, HBCTrojWym, PDwuWym);

    for (int i = 0; i < siatka.nE; ++i) {
        //showTab4x4Dynamic(HBCTrojWym[i]);
    }
    HBCdoRownania = agregacjaHBC(HBCTrojWym, siatka);

    //cout << endl << "Agregacja macierzy HBC:" << endl;
    //for (int i = 0; i < siatka.nN; i++)
    //{
    //    for (int j = 0; j < siatka.nN; j++)
    //    {
    //        std::cout << std::setprecision(4) << "[" << HBCdoRownania[i][j] << "]";
    //    }cout << endl;
    //}



    for (int i = 0; i < siatka.nN; i++)
    {
        for (int j = 0; j < siatka.nN; j++)
        {
            tabGlobal1[i][j] = HBCdoRownania[i][j] + HdoRownania[i][j];
        }
    }

    //cout << "Zagregowana macierz HBC + H w rownaniu jako parametr ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,," << endl;
    //cout << endl << "Agregacja - macierz globalna:" << endl;
    //for (int i = 0; i < siatka.nN; i++) {
    //    for (int j = 0; j < siatka.nN; j++) {
    //        std::cout << std::setprecision(4) << "[" << tabGlobal1[i][j] << "]";
    //    }cout << endl;
    //}


    double* Pglobal = new double[siatka.nN];
    double* PglobalDoRownania = new double[siatka.nN];
    for (int i = 0; i < siatka.nN; i++)
    {
        Pglobal[i] = 0;
        PglobalDoRownania[i] = 0;
    }

    PglobalDoRownania = generycznoscP(siatka, Pglobal, PDwuWym);

    //for (int i = 0; i < siatka.nN; i++) {
    //    for (int j = 0; j < siatka.nN; j++) {
    //        std::cout << std::setprecision(4) << "[" << tabGlobal3[i][j] << "]";
    //    }cout << endl;
    //}

    //cout << "Agregacja matrixH" << endl;
    agregacja(tabGlobal1, siatka);

    //for (int i = 0; i < siatka.nN; i++) {
    //    for (int j = 0; j < siatka.nN; j++) {
    //        std::cout << std::setprecision(4) << "[" << HdoRownania[i][j] << "]";
    //    }cout << endl;
    //}


    double* ResultTabGaus = new double[siatka.nN];
    for (int i = 0; i < siatka.nN; i++)
    {
        ResultTabGaus[i] = 1;
    }

    double* vectorTemp = new double[siatka.nN];
    for (int i = 0; i < siatka.nN; i++) {
        vectorTemp[i] = globalData.initialTemp;
    }


    double l_iter = 2;
    for (int iter = 0; iter < l_iter; iter++)
    {
        for (int n = 0; n < siatka.nN; n++) {
            siatka.tabNodes[n].temp0 = vectorTemp[n];
        }
        cout << endl << "###########################################################################################################" << endl;
        cout << endl << "                                Wyniki dla iteracji numer: " << iter << endl;
        cout << endl << "###########################################################################################################" << endl;
        cout << "Macierz C:" << endl;
        showTab16x16Dynamic(CdoRownania);
        cout << "Macierz H:" << endl;
        showTab16x16Dynamic(HdoRownania);
        cout << "Macierz HBC:" << endl;
        showTab16x16Dynamic(HBCdoRownania);

        calculateMatrixH(globalData, siatka, CdoRownania, HdoRownania, matrixH, tabGlobalH, HBCdoRownania);
        cout << endl << "Global H (H + HBC):" << endl;
        showTab16x16Dynamic(tabGlobalH);
        P_vector(globalData, siatka, CdoRownania, PglobalDoRownania, vectorP);
        calculateHP(siatka, HPTab, matrixH, vectorP);

        vectorTemp = EliminacjaGaussa(HPTab, siatka.nN);

        for (int i = 0; i < siatka.nN; i++) {
            for (int j = 0; j < siatka.nN; j++) {
                HPTab[i][j] = 0;
                matrixH[i][j] = 0;
            }
        }
    }
}


