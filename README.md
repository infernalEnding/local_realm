# local_realm
Societal simulations.

iteration 1:
g++ -o economy_simulation economic_simulation.cpp
./economy_simulation


iteration 2:
g++ -o economy_visualization economy_simulation.cpp -lsfml-graphics -lsfml-window -lsfml-system -lm
./economy_visualization\