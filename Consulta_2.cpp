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
#include <random>


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
    return !(nodo->maxLat < latMin || nodo->minLat > latMax ||
             nodo->maxLon < lonMin || nodo->minLon > lonMax);
}

void buscarEnZona(NodoArbolBola *nodo, double latMin, double latMax, double lonMin, double lonMax, std::vector<int> &resultado)
{
    if (!nodo || !intersectaGeograficamente(nodo, latMin, latMax, lonMin, lonMax))
        return;

    if (!nodo->izquierdo && !nodo->derecho) 
    {
        for (const auto &idx : nodo->indices)
        {
            const Punto &p = puntos[idx];
            if (p.latitud >= latMin && p.latitud <= latMax &&
                p.longitud >= lonMin && p.longitud <= lonMax)
            {
                resultado.push_back(idx);
            }
        }
    }
    else
    {
        buscarEnZona(nodo->izquierdo, latMin, latMax, lonMin, lonMax, resultado);
        buscarEnZona(nodo->derecho, latMin, latMax, lonMin, lonMax, resultado);
    }
}

float distanciaEuclidiana(const std::vector<float> &a, const std::vector<float> &b)
{
    float dist = 0.0f;
    for (size_t i = 0; i < a.size(); ++i)
        dist += (a[i] - b[i]) * (a[i] - b[i]);
    return std::sqrt(dist);
}

void inicializarCentroides(std::vector<std::vector<float>> &centroides, const std::vector<int> &indices, int K, int dim) {
    std::mt19937 gen(123);
    std::uniform_int_distribution<> dis(0, indices.size() - 1);
    for (int k = 0; k < K; ++k) {
        centroides[k] = puntos[indices[dis(gen)]].atributos;
    }
}

void MiniBatchKMeans(const std::vector<int> &indices, int K, int max_iter, int batch_size, std::vector<std::vector<float>> &centroides, std::vector<int> &labels) {
    int N = indices.size();
    int D = puntos[indices[0]].atributos.size();
    centroides.resize(K, std::vector<float>(D));
    labels.resize(N, -1);

    inicializarCentroides(centroides, indices, K, D);
    std::vector<int> counts(K, 0);

    std::mt19937 gen(456);
    std::uniform_int_distribution<> dis(0, N - 1);

    for (int iter = 0; iter < max_iter; ++iter) {
        std::vector<int> batch_idxs(batch_size);
        for (int i = 0; i < batch_size; ++i)
            batch_idxs[i] = dis(gen);

        for (int i = 0; i < batch_size; ++i) {
            int idx_global = indices[batch_idxs[i]];
            const std::vector<float> &attrs = puntos[idx_global].atributos;

            float min_dist = std::numeric_limits<float>::max();
            int best_k = -1;
            for (int k = 0; k < K; ++k) {
                float dist = distanciaEuclidiana(attrs, centroides[k]);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_k = k;
                }
            }

            counts[best_k]++;
            float eta = 1.0f / counts[best_k];
            for (int d = 0; d < D; ++d)
                centroides[best_k][d] = (1 - eta) * centroides[best_k][d] + eta * attrs[d];
        }
    }

    for (int i = 0; i < N; ++i) {
        const std::vector<float> &attrs = puntos[indices[i]].atributos;
        float min_dist = std::numeric_limits<float>::max();
        int best_k = -1;
        for (int k = 0; k < K; ++k) {
            float dist = distanciaEuclidiana(attrs, centroides[k]);
            if (dist < min_dist) {
                min_dist = dist;
                best_k = k;
            }
        }
        labels[i] = best_k;
    }
}

void sugerirParametros(int n, int &K, int &batch_size, int &max_iter) {
    if (n <= 1000) {
        K = std::max(2, static_cast<int>(std::sqrt(n)));
        batch_size = std::max(10, n / 5);
        max_iter = 50;
    } else if (n <= 10000) {
        K = std::sqrt(n / 2);
        batch_size = 200;
        max_iter = 80;
    } else {
        K = 300;
        batch_size = 500;
        max_iter = 100;
    } 
}

std::vector<int> sampleIndices(const std::vector<int>& indices, int sample_size) {
    std::vector<int> muestra = indices;
    std::mt19937 gen(123);
    std::shuffle(muestra.begin(), muestra.end(), gen);
    if (sample_size < muestra.size())
        muestra.resize(sample_size);
    return muestra;
}

int elegirTamañoMuestra(int n) {
    if (n <= 5000) return n;
    if (n <= 100000) return 10000;
    if (n <= 1000000) return 20000;
    return 30000;
}



int main()
{
    std::string archivo_puntos = "puntos_10m.bin";
    std::string archivo_arbol = "arbol_bola_preprocesado10m.bin";

    cargarVectorPuntos(archivo_puntos);
    NodoArbolBola *raiz = cargarArbol(archivo_arbol);

    if (!raiz) return 1;

    std::cout << "Árbol cargado correctamente.\n";
    std::cout << "Puntos cargados: " << puntos.size() << "\n";

    double latMin = -100, latMax = 100;
    double lonMin = -100.0, lonMax = 100.0;

    std::vector<int> resultados;
    auto inicioCargaDatos = std::chrono::high_resolution_clock::now();

    buscarEnZona(raiz, latMin, latMax, lonMin, lonMax, resultados);

    std::cout << "Puntos filtrados: " << resultados.size() << "\n";

    if (!resultados.empty()) {
        int sample_size = elegirTamañoMuestra(resultados.size());
        std::vector<int> muestra = sampleIndices(resultados, sample_size);

        int K, batch_size, max_iter;
        sugerirParametros(sample_size, K, batch_size, max_iter);

        std::cout << "K = " << K << ", batch_size = " << batch_size << ", max_iter = " << max_iter << "\n";

        std::vector<std::vector<float>> centroides;
        std::vector<int> labels_muestra;

        MiniBatchKMeans(muestra, K, max_iter, batch_size, centroides, labels_muestra);

        std::vector<int> labels(resultados.size());
        for (size_t i = 0; i < resultados.size(); ++i) {
            const auto& attrs = puntos[resultados[i]].atributos;
            float min_dist = std::numeric_limits<float>::max();
            int best_k = -1;
            for (int k = 0; k < K; ++k) {
                float dist = distanciaEuclidiana(attrs, centroides[k]);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_k = k;
                }
            }
            labels[i] = best_k;
        }

        auto finCargaDatos = std::chrono::high_resolution_clock::now();

        std::cout << "Tiempo de clusterización: "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(finCargaDatos - inicioCargaDatos).count()
                  << " ms\n";

        std::vector<std::pair<int, int>> puntosOrdenados;
        for (size_t i = 0; i < resultados.size(); ++i) {
            puntosOrdenados.emplace_back(labels[i], resultados[i]);
        }
        std::sort(puntosOrdenados.begin(), puntosOrdenados.end());

        std::unordered_map<int, int> conteo_clusters;
        std::ofstream archivo_salida("clusters_filtrados.csv");
        archivo_salida << "id,lat,lon,cluster\n";

        for (const auto& par : puntosOrdenados) {
            int cluster_id = par.first;
            int idx = par.second;
            const Punto& p = puntos[idx];
            archivo_salida << p.id << "," << p.latitud << "," << p.longitud << "," << cluster_id << "\n";
            conteo_clusters[cluster_id]++;
        }

        archivo_salida.close();
        std::cout << "Archivo clusters_filtrados.csv generado con éxito.\n";

        std::ofstream resumen("resumen_clusters.txt");
        for (const auto& par : conteo_clusters) {
            resumen << "Cluster " << par.first << ": " << par.second << " puntos\n";
        }
        resumen.close();
        std::cout << "Archivo resumen_clusters.txt generado con éxito.\n";

    } else {
        std::cout << "No hay puntos en la región para clusterizar.\n";
    }

    return 0;
}
