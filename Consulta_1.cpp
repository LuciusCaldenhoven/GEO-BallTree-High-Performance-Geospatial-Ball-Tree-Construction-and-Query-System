#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <unordered_map>
#include <queue>
#include <iomanip>


struct Punto
{
    int id;
    double latitud, longitud;
    std::vector<float> atributos;
    Punto(int i, double lat, double lon, const std::vector<float> &attrs)
        : id(i), latitud(lat), longitud(lon), atributos(attrs) {}
    Punto() {}
};
std::vector<Punto> puntos; 
std::unordered_map<int, int> hashMapTripID; 

struct NodoArbolBola
{
    std::vector<float> centro;
    float radio;
    double minLat, maxLat, minLon, maxLon;
    NodoArbolBola *izquierdo = nullptr;
    NodoArbolBola *derecho = nullptr;
    std::vector<int> indices;
};

struct ResultadoSimilitud
{
    float distancia;
    int idxPunto;
    bool operator<(const ResultadoSimilitud &otro) const
    {
        return distancia < otro.distancia;
    }
};


void cargarVectorPuntos(const std::string &nombre)
{
    std::ifstream archivo(nombre, std::ios::binary);
    if (!archivo)
    {
        std::cerr << "No se pudo abrir archivo de puntos\n";
        exit(1);
    }
    size_t n_puntos;
    archivo.read(reinterpret_cast<char*>(&n_puntos), sizeof(size_t));
    puntos.resize(n_puntos);

    for (size_t i = 0; i < n_puntos; ++i) {
        int id;
        double lat, lon;
        archivo.read(reinterpret_cast<char*>(&id), sizeof(int));
        archivo.read(reinterpret_cast<char*>(&lat), sizeof(double));
        archivo.read(reinterpret_cast<char*>(&lon), sizeof(double));
        size_t n_atrib;
        archivo.read(reinterpret_cast<char*>(&n_atrib), sizeof(size_t));
        std::vector<float> atributos(n_atrib);
        archivo.read(reinterpret_cast<char*>(atributos.data()), n_atrib * sizeof(float));
        puntos[i] = Punto(id, lat, lon, atributos);
        hashMapTripID[id] = i;
    }
    archivo.close();
}


NodoArbolBola *cargarNodo(std::ifstream &archivo)
{
    NodoArbolBola *nodo = new NodoArbolBola();

    size_t dim;
    archivo.read(reinterpret_cast<char *>(&dim), sizeof(size_t));
    nodo->centro.resize(dim);
    archivo.read(reinterpret_cast<char *>(nodo->centro.data()), dim * sizeof(float));

    archivo.read(reinterpret_cast<char *>(&nodo->radio), sizeof(float));
    archivo.read(reinterpret_cast<char *>(&nodo->minLat), sizeof(double));
    archivo.read(reinterpret_cast<char *>(&nodo->maxLat), sizeof(double));
    archivo.read(reinterpret_cast<char *>(&nodo->minLon), sizeof(double));
    archivo.read(reinterpret_cast<char *>(&nodo->maxLon), sizeof(double));

    bool esHoja;
    archivo.read(reinterpret_cast<char *>(&esHoja), sizeof(bool));

    if (esHoja)
    {
        size_t n_indices;
        archivo.read(reinterpret_cast<char *>(&n_indices), sizeof(size_t));
        nodo->indices.resize(n_indices);
        archivo.read(reinterpret_cast<char *>(nodo->indices.data()), n_indices * sizeof(int));
    }
    else
    {
        nodo->izquierdo = cargarNodo(archivo);
        nodo->derecho = cargarNodo(archivo);
    }

    return nodo;
}

NodoArbolBola *cargarArbol(const std::string &archivo_nombre)
{
    std::ifstream archivo(archivo_nombre, std::ios::binary);
    if (!archivo)
    {
        std::cerr << "No se pudo abrir el archivo.\n";
        return nullptr;
    }
    NodoArbolBola *raiz = cargarNodo(archivo);
    archivo.close();
    return raiz;
}


bool intersectaGeograficamente(const NodoArbolBola *nodo, double latMin, double latMax, double lonMin, double lonMax)
{
    return !(nodo->maxLat < latMin || nodo->minLat > latMax || nodo->maxLon < lonMin || nodo->minLon > lonMax);
}


float distanciaEuclidiana(const std::vector<float> &a, const std::vector<float> &b)
{
    float dist = 0.0f;
    for (size_t i = 0; i < a.size(); ++i)
        dist += (a[i] - b[i]) * (a[i] - b[i]);
    return std::sqrt(dist);
}

void buscarSimilaresConFiltroGeografico(NodoArbolBola *nodo, const Punto &consulta, double latMin, double latMax, double lonMin, double lonMax, int k,
    std::priority_queue<ResultadoSimilitud> &mejores)
{
    if (!nodo || !intersectaGeograficamente(nodo, latMin, latMax, lonMin, lonMax))
        return;

    float distCentro = distanciaEuclidiana(nodo->centro, consulta.atributos);
    float cotaInferior = distCentro - nodo->radio;
    if (mejores.size() == k && cotaInferior >= mejores.top().distancia)
        return;

    if (!nodo->izquierdo && !nodo->derecho)
    {
        for (const auto &idx : nodo->indices)
        {
            const Punto &p = puntos[idx];
            if (p.latitud >= latMin && p.latitud <= latMax && p.longitud >= lonMin && p.longitud <= lonMax && p.id != consulta.id)
            {
                float dist = distanciaEuclidiana(p.atributos, consulta.atributos);
                if (mejores.size() < k)
                {
                    mejores.push({dist, idx});
                }   
                else if (dist < mejores.top().distancia)
                {
                    mejores.pop();
                    mejores.push({dist, idx});
                }
            }
        }
    }
    else
    {
        NodoArbolBola *primero = nodo->izquierdo;
        NodoArbolBola *segundo = nodo->derecho;
        if (distanciaEuclidiana(nodo->derecho->centro, consulta.atributos) < distanciaEuclidiana(nodo->izquierdo->centro, consulta.atributos))
            std::swap(primero, segundo);

        buscarSimilaresConFiltroGeografico(primero, consulta, latMin, latMax, lonMin, lonMax, k, mejores);
        buscarSimilaresConFiltroGeografico(segundo, consulta, latMin, latMax, lonMin, lonMax, k, mejores);
    }
}




int main()
{
    std::string archivo_puntos = "puntos_10m.bin";
    std::string archivo_arbol = "arbol_bola_preprocesado10m.bin";

    cargarVectorPuntos(archivo_puntos);
    auto inicioCarga = std::chrono::high_resolution_clock::now();
    NodoArbolBola *raiz = cargarArbol(archivo_arbol);
    auto finCarga = std::chrono::high_resolution_clock::now();

    if (!raiz) return 1;

    std::cout << "Árbol cargado correctamente.\n";
    std::cout << "Puntos cargados: " << puntos.size() << "\n";

    int tripID_consulta = 1000000;

    int idxConsulta = hashMapTripID[tripID_consulta];
    const Punto &consulta = puntos[idxConsulta];

    int k = 10000;
    double latMin = -100;
    double latMax = 100;
    double lonMin = -100.0;
    double lonMax = 100;

    std::priority_queue<ResultadoSimilitud> mejores;

    auto inicioCargaDatos = std::chrono::high_resolution_clock::now();
    buscarSimilaresConFiltroGeografico(raiz, consulta, latMin, latMax, lonMin, lonMax, k, mejores);
    auto finCargaDatos = std::chrono::high_resolution_clock::now();
    
    std::cout << "Tiempo: "<< std::chrono::duration_cast<std::chrono::milliseconds>(finCargaDatos - inicioCargaDatos).count()<< " ms\n";


    std::vector<ResultadoSimilitud> ordenados;
    while (!mejores.empty())
    {
        ordenados.push_back(mejores.top());
        mejores.pop();
    }
    std::reverse(ordenados.begin(), ordenados.end());


    std::ofstream archivo("resultados_similares.csv");

    archivo << "Ranking,TripID,Distancia,Latitud,Longitud\n";
    int rank = 1;
    for (const auto &r : ordenados)
    {
        const Punto &p = puntos[r.idxPunto];
        archivo << rank << "," << p.id << "," << std::fixed << std::setprecision(6)<< r.distancia << "," << p.latitud << "," << p.longitud << "\n";
        rank++;
    }
    archivo.close();
    std::cout << "✅ Resultados guardados en resultados_similares.csv\n";

    

    return 0;
}
