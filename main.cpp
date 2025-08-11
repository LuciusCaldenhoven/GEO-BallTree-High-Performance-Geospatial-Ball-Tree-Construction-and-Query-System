#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <unordered_map>


struct Punto
{
    int id;
    double latitud, longitud;
    std::vector<float> atributos;

    Punto(int i, double lat, double lon, const std::vector<float> &attrs)
        : id(i), latitud(lat), longitud(lon), atributos(attrs) {}

    Punto() {}
};

struct NodoArbolBola
{
    std::vector<float> centro;
    float radio;
    double minLat, maxLat, minLon, maxLon;
    NodoArbolBola *izquierdo;
    NodoArbolBola *derecho;
    std::vector<int> indices; 

    NodoArbolBola() : izquierdo(nullptr), derecho(nullptr) {}
};

std::vector<float> calcularCentro(const std::vector<int> &indices, const std::vector<Punto> &puntos)
{
    size_t dim = puntos[indices[0]].atributos.size();
    std::vector<float> centro(dim, 0.0f);
    for (int idx : indices)
        for (size_t i = 0; i < dim; ++i)
            centro[i] += puntos[idx].atributos[i];
    for (float &val : centro)
        val /= indices.size();
    return centro;
}

float distanciaEuclidiana(const std::vector<float> &a, const std::vector<float> &b)
{
    float distancia = 0.0f;
    for (size_t i = 0; i < a.size(); ++i)
        distancia += (a[i] - b[i]) * (a[i] - b[i]);
    return std::sqrt(distancia);
}

float calcularRadio(const std::vector<int> &indices, const std::vector<Punto> &puntos, const std::vector<float> &centro)
{
    float maxDistancia = 0.0f;
    for (int idx : indices)
        maxDistancia = std::max(maxDistancia, distanciaEuclidiana(puntos[idx].atributos, centro));
    return maxDistancia;
}

void calcularBoundingBox(const std::vector<int> &indices, const std::vector<Punto> &puntos, double &minLat, double &maxLat, double &minLon, double &maxLon)
{
    minLat = std::numeric_limits<double>::max();
    maxLat = std::numeric_limits<double>::lowest();
    minLon = std::numeric_limits<double>::max();
    maxLon = std::numeric_limits<double>::lowest();
    for (int idx : indices)
    {
        const Punto &p = puntos[idx];
        minLat = std::min(minLat, p.latitud);
        maxLat = std::max(maxLat, p.latitud);
        minLon = std::min(minLon, p.longitud);
        maxLon = std::max(maxLon, p.longitud);
    }
}


NodoArbolBola *construirArbolBola(const std::vector<Punto> &puntos, std::vector<int> indices, int tamaño_hoja = 256)
{
    NodoArbolBola *nodo = new NodoArbolBola();

    calcularBoundingBox(indices, puntos, nodo->minLat, nodo->maxLat, nodo->minLon, nodo->maxLon);
    nodo->centro = calcularCentro(indices, puntos);
    nodo->radio = calcularRadio(indices, puntos, nodo->centro);

    if (indices.size() <= tamaño_hoja)
    {
        nodo->indices = indices;
        return nodo;
    }
    size_t dimensiones = puntos[indices[0]].atributos.size();
    size_t mejor_dim = 0;
    float mejor_var = -1.0f;
    for (size_t d = 0; d < dimensiones; ++d)
    {
        float media = 0.0f;
        for (int idx : indices)
            media += puntos[idx].atributos[d];
        media /= indices.size();

        float varianza = 0.0f;
        for (int idx : indices)
            varianza += (puntos[idx].atributos[d] - media) * (puntos[idx].atributos[d] - media);
        if (varianza > mejor_var)
        {
            mejor_var = varianza;
            mejor_dim = d;
        }
    }

    std::nth_element(indices.begin(), indices.begin() + indices.size() / 2, indices.end(),[&puntos, mejor_dim](int a, int b)
    {
        return puntos[a].atributos[mejor_dim] < puntos[b].atributos[mejor_dim];
    });

    std::vector<int> izquierda(indices.begin(), indices.begin() + indices.size() / 2);
    std::vector<int> derecha(indices.begin() + indices.size() / 2, indices.end());

    nodo->izquierdo = construirArbolBola(puntos, izquierda, tamaño_hoja);
    nodo->derecho = construirArbolBola(puntos, derecha, tamaño_hoja);

    return nodo;
}


std::vector<Punto> leerCSV(const std::string &archivo)
{
    std::vector<Punto> puntos;
    std::ifstream file(archivo);
    std::string linea;
    std::getline(file, linea); 

    while (std::getline(file, linea))
    {
        std::stringstream ss(linea);
        std::string celda;
        std::vector<float> atributos;
        double lat, lon;
        int tripID = 0;
        int col = 0;
        while (std::getline(ss, celda, ','))
        {
            float val = std::stof(celda);
            if (col == 0)
                tripID = static_cast<int>(val);
            else if (col == 1)
                lat = val;
            else if (col == 2)
                lon = val;
            else if (col >= 3)
                atributos.push_back(val);
            col++;
        }
        puntos.emplace_back(tripID, lat, lon, atributos);
    }
    return puntos;
}


void guardarVectorPuntos(const std::vector<Punto> &puntos, const std::string &nombre) {
    std::ofstream archivo(nombre, std::ios::binary);
    if (!archivo) {
        std::cerr << "No se pudo abrir archivo de puntos para guardar\n";
        return;
    }
    size_t n_puntos = puntos.size();
    archivo.write(reinterpret_cast<const char*>(&n_puntos), sizeof(size_t));
    for (const auto &p : puntos) {
        archivo.write(reinterpret_cast<const char*>(&p.id), sizeof(int));
        archivo.write(reinterpret_cast<const char*>(&p.latitud), sizeof(double));
        archivo.write(reinterpret_cast<const char*>(&p.longitud), sizeof(double));
        size_t n_atrib = p.atributos.size();
        archivo.write(reinterpret_cast<const char*>(&n_atrib), sizeof(size_t));
        archivo.write(reinterpret_cast<const char*>(p.atributos.data()), n_atrib * sizeof(float));
    }
    archivo.close();
    std::cout << "Vector de puntos guardado en: " << nombre << std::endl;
}


void serializarNodo(NodoArbolBola *nodo, std::ofstream &archivo)
{
    size_t dim = nodo->centro.size();
    archivo.write(reinterpret_cast<const char *>(&dim), sizeof(size_t));
    archivo.write(reinterpret_cast<const char *>(nodo->centro.data()), dim * sizeof(float));
    archivo.write(reinterpret_cast<const char *>(&nodo->radio), sizeof(float));
    archivo.write(reinterpret_cast<const char *>(&nodo->minLat), sizeof(double));
    archivo.write(reinterpret_cast<const char *>(&nodo->maxLat), sizeof(double));
    archivo.write(reinterpret_cast<const char *>(&nodo->minLon), sizeof(double));
    archivo.write(reinterpret_cast<const char *>(&nodo->maxLon), sizeof(double));

    bool esHoja = (nodo->izquierdo == nullptr && nodo->derecho == nullptr);
    archivo.write(reinterpret_cast<const char *>(&esHoja), sizeof(bool));

    if (esHoja)
    {
        size_t n_indices = nodo->indices.size();
        archivo.write(reinterpret_cast<const char *>(&n_indices), sizeof(size_t));
        archivo.write(reinterpret_cast<const char *>(nodo->indices.data()), n_indices * sizeof(int));
    }
    else
    {
        serializarNodo(nodo->izquierdo, archivo);
        serializarNodo(nodo->derecho, archivo);
    }
}

void guardarArbol(NodoArbolBola *raiz, const std::string &nombre_archivo)
{
    std::ofstream archivo(nombre_archivo, std::ios::binary);
    if (!archivo)
    {
        std::cerr << "No se pudo abrir el archivo para escritura.\n";
        return;
    }
    serializarNodo(raiz, archivo);
    archivo.close();
    std::cout << "Árbol guardado en: " << nombre_archivo << std::endl;
}


int main()
{
    std::string archivo = "processed_data_complete.csv";

    auto inicioCarga = std::chrono::high_resolution_clock::now();
    std::vector<Punto> puntos = leerCSV(archivo);
    auto finCarga = std::chrono::high_resolution_clock::now();
    std::cout << "Puntos cargados: " << puntos.size() << std::endl;

    guardarVectorPuntos(puntos, "puntos_10m.bin");


    std::vector<int> indices(puntos.size());
    for (int i = 0; i < puntos.size(); ++i) indices[i] = i;

    auto inicioConstruccion = std::chrono::high_resolution_clock::now();
    NodoArbolBola *raiz = construirArbolBola(puntos, indices);
    auto finConstruccion = std::chrono::high_resolution_clock::now();
    std::cout << "Árbol Bola Geográfico construido con éxito.\n";

    auto inicioGuardado = std::chrono::high_resolution_clock::now();
    guardarArbol(raiz, "arbol_bola_preprocesado10m.bin");
    auto finGuardado = std::chrono::high_resolution_clock::now();

    auto tiempoCarga = std::chrono::duration_cast<std::chrono::milliseconds>(finCarga - inicioCarga).count();
    auto tiempoConstruccion = std::chrono::duration_cast<std::chrono::milliseconds>(finConstruccion - inicioConstruccion).count();
    auto tiempoGuardado = std::chrono::duration_cast<std::chrono::milliseconds>(finGuardado - inicioGuardado).count();

    std::cout << "Tiempo de carga: " << tiempoCarga << " ms\n";
    std::cout << "Tiempo de construcción: " << tiempoConstruccion << " ms\n";
    std::cout << "Tiempo de guardado: " << tiempoGuardado << " ms\n";

    return 0;
}


