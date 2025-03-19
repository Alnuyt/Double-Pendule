# Simulation du Double Pendule (Méthode des Différences Finies)
Ce projet, réalisé en groupe dans le cadre du cours **LPHYS1303**, implémente une simulation numérique du **double pendule**, un système dynamique non linéaire connu pour son **comportement chaotique**. L'intégration des équations du mouvement est réalisée avec la **méthode Runge-Kutta d'ordre 4 (RK4)**.
## Description
Le double pendule est constitué de deux masses suspendues l'une à l'autre. Ses équations du mouvement sont obtenues à l'aide du **formalisme lagrangien** et présentent un **couplage non linéaire** rendant le système très sensible aux conditions initiales.

L'objectif est :
- De résoudre numériquement ces équations avec **Runge-Kutta 4**
- D'analyser qualitativement différentes trajectoires
- D'étudier la **sensibilité aux conditions initiales** et le **chaos** à travers l'exposant de Lyapunov
## Utilisation 
Exécutez le code pour lancer la simulation : [Pendule_Double.py](Pendule_Double.py)
## Résultats
- **Simulation des trajectoires :** analyse du mouvement des masses pour différentes conditions initiales.
- **Mise en évidence du chaos :** visualisation de l'évolution des angles et vitesses.
- **Analyse quantitative :** étude de l'exposant de **Lyapunov** pour caractériser la divergence des trajectoires.

