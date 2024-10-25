#include <iostream>
#include <vector>
#include <cmath>
#include <string>
using namespace std;

class Node {
public:
    string name;
    double A;      // Total factor productivity
    double K;      // Capital stock
    double N;      // Population
    double s;      // Savings rate
    double alpha;  // Output elasticity of capital
    double delta;  // Depreciation rate
    double Y;      // Output
    double C;      // Consumption
    double S;      // Savings
    double I;      // Investment
    double NX;     // Net exports

    Node(string n, double A_init, double K_init, double N_init, double s_init,
         double alpha_init, double delta_init)
        : name(n), A(A_init), K(K_init), N(N_init), s(s_init),
          alpha(alpha_init), delta(delta_init), Y(0), C(0),
          S(0), I(0), NX(0) {}

    // Production function
    void produce() {
        Y = A * pow(K, alpha) * pow(N, 1 - alpha);
    }

    // Update capital stock
    void updateCapital() {
        K = (1 - delta) * K + I;
    }

    // Update population
    void updatePopulation(double n_growth) {
        N = N + n_growth * N;
    }

    // Compute savings and investment
    void computeSavingsAndInvestment() {
        S = s * Y;
        I = S + NX;
    }

    // Compute consumption
    void consume() {
        C = (1 - s) * Y;
    }
};



double computeTradeFlow(Node &node_i, Node &node_j, double distance,
                        double G, double beta) {
    return G * (node_i.Y * node_j.Y) / pow(distance, beta);
}




int main() {
    // Simulation parameters
    const int TIME_STEPS = 10;
    const double G = 1.0;     // Gravity model constant
    const double beta = 1.0;  // Distance decay parameter
    const double n_0 = 0.01;  // Base population growth rate
    const double eta = 0.0001; // Sensitivity of population growth
    const double y_star = 1.0; // Subsistence income level

    // Initialize nodes
    Node A("A", 1.0, 100.0, 1000.0, 0.2, 0.3, 0.05);
    Node B("B", 1.0, 80.0, 800.0, 0.25, 0.3, 0.05);
    Node C("C", 1.0, 70.0, 700.0, 0.22, 0.3, 0.05);
    Node D("D", 1.0, 60.0, 600.0, 0.18, 0.3, 0.05);
    Node E("E", 1.0, 50.0, 500.0, 0.20, 0.3, 0.05);

    vector<Node> nodes = {A, B, C, D, E};

    // Distance matrix (symmetric)
    vector<vector<double>> distances = {
        {0, 10, 20, 30, 40},
        {10, 0, 15, 25, 35},
        {20, 15, 0, 12, 22},
        {30, 25, 12, 0, 14},
        {40, 35, 22, 14, 0}
    };

    // Simulation loop
    for (int t = 0; t < TIME_STEPS; ++t) {
        // Step 1: Production
        for (auto &node : nodes) {
            node.produce();
        }

        // Step 2: Trade between nodes
        vector<vector<double>> trade_flows(nodes.size(),
                                           vector<double>(nodes.size(), 0));
        for (size_t i = 0; i < nodes.size(); ++i) {
            for (size_t j = 0; j < nodes.size(); ++j) {
                if (i != j) {
                    trade_flows[i][j] = computeTradeFlow(
                        nodes[i], nodes[j], distances[i][j], G, beta);
                }
            }
        }

        // Step 3: Net exports
        for (size_t i = 0; i < nodes.size(); ++i) {
            double total_exports = 0.0;
            double total_imports = 0.0;
            for (size_t j = 0; j < nodes.size(); ++j) {
                if (i != j) {
                    total_exports += trade_flows[i][j];
                    total_imports += trade_flows[j][i];
                }
            }
            nodes[i].NX = total_exports - total_imports;
        }

        // Step 4: Consumption and Savings
        for (auto &node : nodes) {
            node.computeSavingsAndInvestment();
            node.consume();
        }

        // Step 5: Capital accumulation
        for (auto &node : nodes) {
            node.updateCapital();
        }

        // Step 6: Population growth
        for (auto &node : nodes) {
            double per_capita_income = node.Y / node.N;
            double n_growth = n_0 + eta * (per_capita_income - y_star);
            node.updatePopulation(n_growth);
        }

        // Output results
        cout << "Time Step " << t + 1 << ":\n";
        for (auto &node : nodes) {
            cout << "Node " << node.name
                 << " - Y: " << node.Y
                 << ", K: " << node.K
                 << ", N: " << node.N
                 << ", NX: " << node.NX << endl;
        }
        cout << "----------------------\n";
    }

    return 0;
}

