# MOMeta_M2_ORO
Projet de Métaheuristique Multi-Objectif M2 ORO

## Comment exécuter ?
```bash
julia main.jl
```
ou à partir du REPL Julia :
```julia
include("main.jl")
```

## Scatter search
Exemple d'exécution du scatter search sur des instances se trouvant dans
un répertoire `instances`, en enregistrant les résultats dans un dossier
`resultats`, avec une population initiale de 200 solutions et le paramètre
`α=0.7`:
```julia
ScatterSearch("instances", "resultats", 200, α=0.7)
```

Les fonctions `runScatterSearch` et `runScatterSearchBattery` sont des
fonctions additionnelles qui permettent de lancer le scatter search plusieurs
fois sur les instances et de sauvegarder les résultats (moyennes) dans un
fichier CSV, et respectivement, de lancer le scatter search sur plusieurs
valeurs de `α` et `TL` et de sauvegarder les résultats dans des fichiers CSV.

E.g:
```julia
runScatterSearch(5, "instances", "resultats", 200, α=0.7)
runScatterSearchBattery(4, 5, "instances", "resultats", 200, αmin=0.5, αmax=0.9, TLmin=3, TLmax=14)
```
Ces appels de fonction vont lancer le scatter search 5 fois sur chaque instance
et sauvegarder les résultats dans un fichier CSV. De plus, la fonction `runScatterSearchBattery`
va lancer le scatter search sur 4 valeurs de `α` et 4 valeurs de `TL` (soit 16 combinaisons).

L'ensemble des paramètres du Scatter Search sont:
- `path`: chemin vers les instances
- `savepath`: chemin pour sauvegarder les résultats
- `popSize`: taille de la population initiale
- `α`: paramètre pour le GRASP, permet de contrôler le caractère aléatoire des
solutions gloutonnes (valeur entre 0 et 1, par défaut 0.7)
- `p`: paramètre pour le GRASP, proportion de concentrateurs de niveau 1 qui
  seront sélectionnés dans la solution gloutonne (valeur entre 0 et 1, par défaut 0.4)
- `TL`: paramètre pour la recherche tabou, taille de la liste tabou
- `kp`: paramètre pour la recherche tabou, permet de régler `k` le nombre
  maximum d'itérations sans amélioration avant d'arrêter la recherche tabou,
  ici `k = kp*nb_concentrateur_niv1` (valeur entre 0 et 1, par défaut 0.4)
- `verboseLevel`: niveau de verbosité, 0 pour aucune verbosité, 1 pour verbosité simple,
  2 pour verbosité complète (incluant les messages de verbosité du solveur utilisé)
- `exact`: paramètre pour le modèle exact, booléen, true pour résoudre le modèle exact
- `unique`: paramètre pour le modèle exact, booléen, true pour arrêter la
  résolution dés qu'une solution non-dominée est trouvée
- `timeout`: paramètre pour le modèle exact, timeout pour la résolution du modèle exact
