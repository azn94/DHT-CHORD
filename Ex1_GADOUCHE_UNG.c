/* GADOUCHE LOTFI 3673119 */
/* UNG RICHARD    3680881 */

#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

/* Valeur à entrer */

#define M               6  // Nombre de bits sur lequel toute clé est encodée, par défaut on a pris 6 comme dans l'exemple du TD
#define NB_SITE         10 // Nombre de pairs, par défaut on a pris 10 comme dans l'exemple du TD

/* TAG des messages */

#define TAG_ID_CHORD    0
#define TAG_FINGER      1
#define TAG_RANK        2

#define TAG_SPECIAL     3
#define TAG_END         4


/* Variable */

int id;                         // id CHORD
int finger_p[M];                // tableau des M fingers du processus
int finger_MPI_rank_p[M];       // tableau des M rangs MPI responsables des M fingers
int succ_p;                     // l’id CHORD de son successeur dans l’anneau
int recherche[2] = {-1,-1};     // recherche[0] := cle_recherché; recherche[1] := id chord responsable de la clé recherché

/* Index recherche */

#define CLE_RECHERCHE   0
#define RESPONSABLE_CLE 1

/* Fonctions */

/*
Fonction pow de math fonctionnait pas donc on a fait nôtre propre fonction power
*/

int power(int x, int y){
    int tmp;
    if( y == 0)
        return 1;
    tmp = power(x, y/2);
    if (y%2 == 0)
        return tmp*tmp;
    else
        return x*tmp*tmp;
}

/*
Fonction pour comparer des éléments deux à deux.
*/
int fonction_cmp ( const void * premier, const void * second ) {
    int premierInt = * (const int *) premier;
    int secondInt = * (const int *) second;
    if(premierInt > secondInt)  return 1;
    if(premierInt == secondInt)  return 0;
    return -1;
}

/*
Fonction qui tire NB_SITE valeurs dans l'ensemble {0, ..., 2^(M) - 1}
pour constituer l'ensemble PI
*/
void AleatoireIDChord(int * tableau_id){
    int max_valeur = power(2, M);           // Valeur maximum qu'un identifiant peut prendre
    int est_id[max_valeur];                 // tableau de booléen qui indique si un id a déjà été choisit comme id CHORD

    for(int i = 0; i < max_valeur; ++i) est_id[i] = 0;

    for (int cpt = 0; cpt < NB_SITE; ++cpt) {
        int id_chord;                       // valeur des id CHORD que l'on va tirer aléatoirement

        // Tant qu'une valeur d'identifiant tiré aléatoirement a déjà été choisit, on réitère
        do{
            id_chord = rand() % max_valeur;
        }while(est_id[id_chord]);

        // Quant on arrive ici, cela signifie qu'un id non choisit a été sélectionné, on l'ajoute au tableau
        est_id[id_chord] = 1;
        tableau_id[cpt] = id_chord;
    }

    // D'après l'énoncé:
    // L’ensemble Π des pairs est arrangé en anneau dans l’ordre croissant  
    // de leur identifiant (dans le sens des aiguilles d’une montre).
    qsort(tableau_id, NB_SITE, sizeof(int),fonction_cmp);
}

/*
Fonction booléenne qui vérifie que k appartient à l'intervalle [a,b[

[a,b[ = {a, a+1, ..., b-1} si a < b
      = {a, a+1, ..., max - 1, 0, 1, ..., b-1} si a >= b
*/
int app(int k, int a, int b, int max) {
if (a < b) return ((k >= a) && (k < b));
return (((k >= 0) && (k < b)) || ((k >= a) && (k < max)));
}

/*
Le processus simulateur se contentera de tirer aléatoirement 
les identifiants CHORD des pairs puis de construire un anneau 
dans lequel les rang MPI seront ordonnés et de déterminer
aléatoirement un ensemble non vide d’initiateurs parmi ces pairs.
*/
void simulateur(int rang, MPI_Status status){

    /*
    Question 1 : Initialisation de la DHT
    */

    int tab_id_chord[NB_SITE];              // Tableau contenant les identités CHORD des NB_SITE pairs, PI.
    //int tab_id_chord[NB_SITE] = {2, 7, 13, 14, 21, 38, 42, 48, 51, 59}; // TEST TD -> mettre en commentaire la ligne précédente
    int deuxi[M];                           // Tableau contenant les valeurs de 2^i
    int max_valeur = power(2, M);           // Valeur maximum qu'un identifiant peut prendre

    
    // On tire aléatoirement les identifiants CHORD des pairs
    AleatoireIDChord(tab_id_chord);         // TEST TD -> mettre en commentaire cette ligne
    
    // On remplit le tableau
    for(int i = 0; i < M; ++i){
        deuxi[i] = power(2,i);
    }

    // Pour chaque processus 
    for(int cpt_proc = 0; cpt_proc < NB_SITE; ++cpt_proc) {
        int id_chord = tab_id_chord[cpt_proc];  // Valeur de l'id CHORD du cpt_proc^ième processus
        int finger[M];                          // tableau contenant les M fingers du cpt_proc^ième processus
        int finger_MPI_rank[M];                 // tableau contenant les M rangs MPI responsables des M fingers

        // Pour chaque finger
        for (int i = 0; i < M; ++i) {
            int valeur = (id_chord + deuxi[i]) % max_valeur;     // (id_chord + 2^i) mod 2^M
            int cpt_id = 0;                                      // Compteur du tableau des id CHORD

            
            // On cherche le plus petit id_chord plus grand que valeur
            while(cpt_id < NB_SITE && valeur > tab_id_chord[cpt_id]){
                ++cpt_id;
            }

            if(cpt_id == NB_SITE) cpt_id = 0;                   // Si cpt_id == NB_SITE, alors valeur est compris entre le plus grand id_chord est 0, donc son responsable est le plus petit id_chord

            // On a trouvé le plus petit id_chord plus grand que valeur
            finger[i] = tab_id_chord[cpt_id];
            finger_MPI_rank[i] = cpt_id + 1;
            
        }

        // cpt_proc + 1 pour ne pas l'envoyé au processus 0 (simulateur) 
        // Envoie des id chord, finger table et leur id MPI associés
        MPI_Send(&id_chord, 1, MPI_INT, cpt_proc + 1, TAG_ID_CHORD, MPI_COMM_WORLD);
        MPI_Send(finger, M, MPI_INT, cpt_proc + 1, TAG_FINGER, MPI_COMM_WORLD);
        MPI_Send(finger_MPI_rank, M, MPI_INT, cpt_proc + 1, TAG_RANK, MPI_COMM_WORLD);
    }

    /*
    Question 2 : Recherche d'une clé
    */

    int rang_initiateur = rand() % NB_SITE;
    recherche[CLE_RECHERCHE] = rand() % max_valeur;

    // TEST TD -> retirer les doubles barres des deux lignes suivantes
    //rang_initiateur = 2; // 2 : index de l'id chord 7, 8 : index de l'id chord 51
    //recherche[CLE_RECHERCHE] = 30; // 30, 0, 1O, 50, 22

    printf("P%d : Le pair p%d recherche la clé %d.\n",rang, tab_id_chord[rang_initiateur], recherche[CLE_RECHERCHE]);
    MPI_Send(recherche, 2, MPI_INT, rang_initiateur + 1, TAG_SPECIAL, MPI_COMM_WORLD);
    // Attente de la reponse du prédecesseur du responsable de la clé recherché
    MPI_Recv(recherche, 2, MPI_INT, MPI_ANY_SOURCE, TAG_SPECIAL, MPI_COMM_WORLD, &status);
    printf("P%d : msg reçu de P%d, le pair p%d est responsable de la clé %d.\n", rang, status.MPI_SOURCE, recherche[RESPONSABLE_CLE],recherche[CLE_RECHERCHE]);

    // TEST de la valeur calculé par les processus pairs

    int val=-1; // Correspond a l'id chord responsable de la clé recherché selon le simulateur

    for (int cpt_finger = 0; cpt_finger < NB_SITE; ++cpt_finger) {
        if (tab_id_chord[cpt_finger] >= recherche[CLE_RECHERCHE]) {
           val = tab_id_chord[cpt_finger];
           break;
        }
    }

	if(val==-1) val=tab_id_chord[0];

   if (val == recherche[RESPONSABLE_CLE]) {
       printf("------------------------------------------BON RESULTAT: processus pair %d, et simulateur trouve %d.------------------------------------------\n", recherche[RESPONSABLE_CLE],val);
   }
   else {
       printf("------------------------------------------MAUVAIS RESULTAT : processus pair %d, et simulateur trouve %d.------------------------------------------\n", recherche[RESPONSABLE_CLE],val);
   }

   // Fin du programme : le simulateur envoie les msgs de fin
   for (int i = 0; i < NB_SITE; ++i) {
       printf("P%d : envoie msg END à P%d.\n",rang,i+1);
       MPI_Send(recherche, 2, MPI_INT, i + 1, TAG_END, MPI_COMM_WORLD);
   }


}

void processus_pair(int rang, MPI_Status status){

    int max_valeur = power(2, M);     // Valeur maximum qu'un identifiant peut prendre

    // Reception des msgs
    MPI_Recv(&id, 1, MPI_INT, 0, TAG_ID_CHORD, MPI_COMM_WORLD, &status);
    MPI_Recv(finger_p, M, MPI_INT, 0, TAG_FINGER, MPI_COMM_WORLD, &status);
    MPI_Recv(finger_MPI_rank_p, M, MPI_INT, 0, TAG_RANK, MPI_COMM_WORLD, &status);

    succ_p = finger_p[0];

    printf("P%d : mon id CHORD est %d.\n",rang,id);
    printf("P%d : mon successeur CHORD est %d.\n",rang,succ_p);
    printf("P%d : Finger table avec M = %d.\n",rang, M);
    printf("P%d : i\t id chord supérieur\n", rang);
    for(int cpt = 0; cpt < M; ++cpt){
        printf("P%d : %d\t %d.\n",rang, cpt, finger_p[cpt]);
    }
    
    // Tant que le simulateur n'a pas envoyé le msg de terminaison on continue
    while(status.MPI_TAG != TAG_END){
        MPI_Recv(recherche, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        // Si on reçoit le msg de terminaison du simulateur
        if(status.MPI_TAG == TAG_END){
            printf("P%d : le processus se termine.\n",rang);
            continue;
        }

        int bool_trouve = 0;    // booléen pour savoir si la clé appartient à un intervalle ]finger_p,id_chord] du processus
        int index_trouve;       // index de la plus grande clé de la finger table tel que la clé recherché est comprise dans ]finger_p[index_trouve],id]  

        
        for(int j = M-1; j>=0; --j ){
            // clé recherché appartient à ]finger_p[j],id_chord] 
            // <=> clé recherché appartient à [finger_p[j]+1,id_chord]
            // <=> clé recherché appartient à [finger_p[j]+1,id_chord+1[
            // car app(k,a,b) vérifie que k appartient à [a,b[
            if(app(recherche[CLE_RECHERCHE],finger_p[j]+1, id+1, max_valeur)){
                bool_trouve = 1;
                index_trouve = j;
                break;
            }
        }

        // Réitère la recherche dans finger_p[index_trouve]
        if(bool_trouve){
            printf("P%d : id chord %d transmet la requête de recherche de la clé %d à id chord %d (P%d).\n",
            rang, id, recherche[CLE_RECHERCHE], finger_p[index_trouve], finger_MPI_rank_p[index_trouve]);
           MPI_Send(recherche, 2, MPI_INT, finger_MPI_rank_p[index_trouve], TAG_SPECIAL, MPI_COMM_WORLD);
        }else{ // le successeur est responsable de la clé recherché
            recherche[RESPONSABLE_CLE] = succ_p;
            printf("P%d : le successeur id chord %d (P%d) est responsable de la clé %d.\n",
            rang, succ_p, finger_MPI_rank_p[0], recherche[CLE_RECHERCHE]);
            MPI_Send(recherche, 2, MPI_INT, 0, TAG_SPECIAL, MPI_COMM_WORLD);
        }
        

    }
}

int main (int argc, char* argv[]) {
    int nb_proc,rang;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);

    if (nb_proc != NB_SITE+1) {
        printf("Nombre de processus incorrect !\n");
        MPI_Finalize();
        exit(2);
    }

    srand (time (NULL));

    MPI_Comm_rank(MPI_COMM_WORLD, &rang);
    MPI_Status status;
  
    if (rang == 0) {
        simulateur(rang,status);
    } else {
        processus_pair(rang,status);
    }
  
    MPI_Finalize();
    return 0;
}
