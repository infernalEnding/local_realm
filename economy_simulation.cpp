// economy_simulation.cpp

#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

using namespace std;

// Node class representing each city-state
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
    float x, y;    // Position on the screen

    Node(string n, double A_init, double K_init, double N_init, double s_init,
         double alpha_init, double delta_init, float x_pos, float y_pos)
        : name(n), A(A_init), K(K_init), N(N_init), s(s_init),
          alpha(alpha_init), delta(delta_init), Y(0), C(0),
          S(0), I(0), NX(0), x(x_pos), y(y_pos) {}

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

// Function to compute trade flow between two nodes
double computeTradeFlow(Node &node_i, Node &node_j, double distance,
                        double G, double beta) {
    return G * (node_i.Y * node_j.Y) / pow(distance, beta);
}

// Function to draw nodes
void drawNodes(sf::RenderWindow &window, const vector<Node> &nodes, double maxGDP) {
    for (const auto &node : nodes) {
        // Create a circle shape for the node
        sf::CircleShape circle;

        // Set the radius based on population
        float radius = static_cast<float>(sqrt(node.N) / 10.0f);
        circle.setRadius(radius);

        // Set the position (subtract radius to center the circle)
        circle.setPosition(node.x - radius, node.y - radius);

        // Set the color based on GDP
        float gdpIntensity = static_cast<float>(node.Y) / maxGDP;
        gdpIntensity = std::min(std::max(gdpIntensity, 0.0f), 1.0f); // Clamp between 0 and 1
        circle.setFillColor(sf::Color(255 * gdpIntensity, 100, 255 * (1 - gdpIntensity)));

        // Draw the circle
        window.draw(circle);
    }
}

// Function to draw trade connections
void drawTradeConnections(sf::RenderWindow &window, const vector<Node> &nodes, const vector<vector<double>> &trade_flows, double maxTradeVolume) {
    for (size_t i = 0; i < nodes.size(); ++i) {
        for (size_t j = i + 1; j < nodes.size(); ++j) {
            double tradeVolume = trade_flows[i][j];

            // Skip if trade volume is negligible
            if (tradeVolume < 1e-2) continue;

            // Create a line between nodes
            sf::Vertex line[] =
            {
                sf::Vertex(sf::Vector2f(nodes[i].x, nodes[i].y)),
                sf::Vertex(sf::Vector2f(nodes[j].x, nodes[j].y))
            };

            // Set the color based on trade volume
            float tradeIntensity = static_cast<float>(tradeVolume) / maxTradeVolume;
            tradeIntensity = std::min(std::max(tradeIntensity, 0.0f), 1.0f); // Clamp between 0 and 1
            sf::Color lineColor(255 * tradeIntensity, 255 * (1 - tradeIntensity), 0);
            line[0].color = lineColor;
            line[1].color = lineColor;

            // Draw the line
            window.draw(line, 2, sf::Lines);
        }
    }
}

// Function to draw node labels
void drawNodeLabels(sf::RenderWindow &window, const vector<Node> &nodes, sf::Font &font) {
    for (const auto &node : nodes) {
        sf::Text label;
        label.setFont(font);
        label.setCharacterSize(12);
        label.setFillColor(sf::Color::White);
        label.setString(node.name);
        label.setPosition(node.x + 5, node.y + 5);

        window.draw(label);
    }
}

int main() {
    // Simulation parameters
    const int TIME_STEPS = 1000; // Adjust as needed
    const double G = 1.0;        // Gravity model constant
    const double beta = 1.0;     // Distance decay parameter
    const double n_0 = 0.01;     // Base population growth rate
    const double eta = 0.0001;   // Sensitivity of population growth
    const double y_star = 1.0;   // Subsistence income level

    // Window dimensions
    const int WINDOW_WIDTH = 800;
    const int WINDOW_HEIGHT = 600;

    // Create the SFML window
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8; // Enable anti-aliasing
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Economic Simulation Visualization", sf::Style::Default, settings);

    // Load a font (make sure you have the font file in your project directory)
    sf::Font font;
    if (!font.loadFromFile("arial.ttf")) {
        std::cerr << "Failed to load font 'arial.ttf'\n";
        // Handle error or exit
    }

    // Initialize nodes with positions
    Node A("A", 1.0, 100.0, 1000.0, 0.2, 0.3, 0.05, 100.0f, 100.0f);
    Node B("B", 1.0, 80.0, 800.0, 0.25, 0.3, 0.05, 300.0f, 100.0f);
    Node C("C", 1.0, 70.0, 700.0, 0.22, 0.3, 0.05, 500.0f, 100.0f);
    Node D("D", 1.0, 60.0, 600.0, 0.18, 0.3, 0.05, 200.0f, 300.0f);
    Node E("E", 1.0, 50.0, 500.0, 0.20, 0.3, 0.05, 400.0f, 300.0f);

    vector<Node> nodes = {A, B, C, D, E};

    // Distance matrix (symmetric)
    vector<vector<double>> distances = {
        {0, 10, 20, 30, 40},
        {10, 0, 15, 25, 35},
        {20, 15, 0, 12, 22},
        {30, 25, 12, 0, 14},
        {40, 35, 22, 14, 0}
    };

    // Initialize trade_flows and maxTradeVolume outside the simulation loop
    vector<vector<double>> trade_flows(nodes.size(), vector<double>(nodes.size(), 0));
    double maxTradeVolume = 0.0;

    // Set a frame rate limit
    window.setFramerateLimit(60);

    // Simulation loop with time control
    sf::Clock clock;
    float simulationTime = 0.0f;
    const float simulationStep = 0.1f; // Adjust the simulation speed

    // Main simulation loop
    while (window.isOpen()) {
        // Handle events
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
            // Handle key presses
            else if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Escape)
                    window.close();
                // Add key controls to interact with the simulation
                else if (event.key.code == sf::Keyboard::Up) {
                    // Increase savings rate for Node A
                    nodes[0].s += 0.01;
                    if (nodes[0].s > 1.0) nodes[0].s = 1.0;
                }
                else if (event.key.code == sf::Keyboard::Down) {
                    // Decrease savings rate for Node A
                    nodes[0].s -= 0.01;
                    if (nodes[0].s < 0.0) nodes[0].s = 0.0;
                }
            }
        }

        // Get elapsed time
        float deltaTime = clock.restart().asSeconds();
        simulationTime += deltaTime;

        // Run simulation step when enough time has passed
        if (simulationTime >= simulationStep) {
            // Step 1: Production
            for (auto &node : nodes) {
                node.produce();
            }

            // Step 2: Trade between nodes
            maxTradeVolume = 0.0; // Reset maxTradeVolume for each simulation step
            for (size_t i = 0; i < nodes.size(); ++i) {
                for (size_t j = 0; j < nodes.size(); ++j) {
                    if (i != j) {
                        trade_flows[i][j] = computeTradeFlow(
                            nodes[i], nodes[j], distances[i][j], G, beta);
                        maxTradeVolume = std::max(maxTradeVolume, trade_flows[i][j]);
                    } else {
                        trade_flows[i][j] = 0.0;
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

            simulationTime = 0.0f; // Reset the timer
        }

        // Compute maximum GDP for scaling
        double maxGDP = 0.0;
        for (const auto &node : nodes) {
            maxGDP = std::max(maxGDP, node.Y);
        }

        // Clear the window
        window.clear(sf::Color::Black);

        // Draw trade connections
        drawTradeConnections(window, nodes, trade_flows, maxTradeVolume);

        // Draw nodes
        drawNodes(window, nodes, maxGDP);

        // Draw node labels
        drawNodeLabels(window, nodes, font);

        // Display the window contents
        window.display();
    }

    return 0;
}