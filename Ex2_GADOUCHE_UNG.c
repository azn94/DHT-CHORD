/* GADOUCHE LOTFI 3673119 */
/* UNG RICHARD    3680881 */


#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

/* Valeur à entrer */

#define M               6   // Nombre de bits sur lequel toute clé est encodée, par défaut on a pris 6 comme dans l'exemple du TD
#define NB_SITE         10  // Nombre de pairs, par défaut on a pris 10 comme dans l'exemple du TD

/* TAG simulateur & création de l'anneau */

#define TAG_ID          0
#define TAG_VOISINS     1
#define TAG_INITIATEUR  2

/* TAG de l'algorithme d'élection d'un leader d'Hirschberg & Sinclair */

#define TAG_IN          3
#define TAG_OUT         4
#define TAG_LEADER      5

/* TAG transmission du tableau des ID CHORDs */

#define TAG_TAB_ID      6
#define TAG_END         7  

/* Variables locales */

int id;                     // id chord du noeud
int initiateur;             // booléen indiquant si le noeud est initiateur dans l'anneau
int voisins[2];             // Tableau qui contient le rang MPI de son prédécesseur (prev) et de son successeur (next) 
int * tab_id;               // Tableau de tous les id chord de l'anneau
int finger[M];              // Tableau des fingers

/* Index voisins */

#define PREV 0
#define NEXT 1

/* Variables Hirschberg & Sinclair */

int initiateur_election;    // booléen indiquant si le noeud est initiateur dans l'élection
int etape   = 0;            // Nombre de tour de l'algorithme d'Hirschberg & Sinclair
int taille  = 0;            // taille de l'anneau que l'on déduira de l'algorithme d'Hirschberg & Sinclair
int cpt_in  = 0;            // Compteur de message in reçu
int cpt_out = 0;            // Compteur de message out reçu
int state   = -1;           // situation vis-à-vis de l'élection (voir les états plus bas)
int msg[2];                 // Structure des messages échangés dans l'algo d'élection

/* Index msg */

#define INIT 0
#define DIST 1

/* STATE*/

#define NSP    -1
#define BATTU   0
#define ELU     1

/* Fonctions simulation et création de l'anneau */

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
}

/*
Fonction qui construit un tableau contenant les prédécesseurs et successeurs de chaque noeud de l'anneau
*/

void ConstruireAnneau(int anneau[][2]){
    for(int rang_mpi = 0; rang_mpi < NB_SITE + 1; ++rang_mpi){
        switch(rang_mpi){
            case 0:
                // On ne fait rien de concret, le noeud O n'est pas dans l'anneau (simulateur)
                anneau[rang_mpi][PREV] = -1;
                anneau[rang_mpi][NEXT] = -1;
                break;

            case 1:
                anneau[rang_mpi][PREV] = NB_SITE;
                anneau[rang_mpi][NEXT] = rang_mpi + 1; // anneau[rang_mpi][NEXT] = 2
                break;

            case NB_SITE:
                anneau[rang_mpi][PREV] = NB_SITE - 1;
                anneau[rang_mpi][NEXT] = 1;
                break;

            default:
                anneau[rang_mpi][PREV] = rang_mpi - 1;
                anneau[rang_mpi][NEXT] = rang_mpi + 1;
        }

    }
}

/*
Fonction qui tire aléatoirement un ensemble non vide d’initiateurs 
parmi ces pairs et retourne le nombre d'initiateur.
*/

int InitiateurAleatoire(int * initiateurs){
    int nbr_initiateur = (rand()%NB_SITE) + 1; // Nombre d'initiateurs non nul

    // Instancie l'ensemble des valeurs du tableau de booléen à faux
    for (int rang_mpi = 0; rang_mpi < NB_SITE + 1; ++rang_mpi) initiateurs[rang_mpi] = 0;

    // À chaque tour de boucle, on tire aléatoirement un nouvel initiateur
    for (int rang_mpi = 0; rang_mpi < nbr_initiateur; ++rang_mpi) {
        int rang_init;  // rang MPI du nouvel initiateur
        
        do{
            // On ajoute + 1 car 0 n'est pas dans l'anneau
            rang_init = (rand()%NB_SITE) + 1;
        }while (initiateurs[rang_init]);

        // En sortant du do_while, on assure que rang_init est un nouvel initiateur
        initiateurs[rang_init] = 1;
    }

    return nbr_initiateur;
}

/*
Le processus simulateur se contentera de tirer aléatoirement 
les identifiants CHORD des pairs puis de construire un anneau 
dans lequel les rang MPI seront ordonnés et de déterminer
aléatoirement un ensemble non vide d’initiateurs parmi ces pairs.
*/
void simulateur(int rang){

    int tab_id_chord[NB_SITE];              // Tableau contenant les identités CHORD des NB_SITE pairs, PI.
    int anneau[NB_SITE+1][2];               // Tableau contenant les prédécesseurs et successeurs de chaque noeud de l'anneau
    int initiateurs[NB_SITE+1];             // Tableau de booléen qui indique si un noeud dans l'anneau est initiateur ou non
    int nbr_initiateur;                     // Nombre d'initiateurs non nul

    // On tire aléatoirement les identifiants CHORD des pairs
    AleatoireIDChord(tab_id_chord);

    // On construit l'anneau dans lequel les rang MPI seront ordonnés
    ConstruireAnneau(anneau);

    // On détermine aléatoirement un ensemble non vide d’initiateurs parmi ces pairs
    nbr_initiateur = InitiateurAleatoire(initiateurs);
    printf("P%d : il y a %d initiateur(s) dans l'anneau.\n", rang, nbr_initiateur);

    // On envoie à chaque noeud de l'anneau son id chord, 
    // son prédécesseur et successeur, ainsi que son état d'initiateur ou non
    for(int rang_mpi = 1; rang_mpi < NB_SITE + 1; ++rang_mpi){
        MPI_Send(&tab_id_chord[rang_mpi - 1], 1, MPI_INT, rang_mpi, TAG_ID, MPI_COMM_WORLD);
        MPI_Send(anneau[rang_mpi], 2, MPI_INT, rang_mpi, TAG_VOISINS, MPI_COMM_WORLD);
        MPI_Send(&initiateurs[rang_mpi], 1, MPI_INT, rang_mpi, TAG_INITIATEUR, MPI_COMM_WORLD);
    }
}

/*
Fonction qui récupère des msgs du simulateur l'id chord, le prédécesseur et 
successeur du noeud, ainsi que son état (si il est initialisateur ou non)
*/

void Instanciation(int rang, MPI_Status status){

    MPI_Recv(&id, 1, MPI_INT, 0, TAG_ID, MPI_COMM_WORLD, &status);
    MPI_Recv(voisins, 2, MPI_INT, 0, TAG_VOISINS, MPI_COMM_WORLD, &status);
    MPI_Recv(&initiateur, 1, MPI_INT, 0, TAG_INITIATEUR, MPI_COMM_WORLD, &status);

    printf("p%d : mon rang MPI est %d.\n", id, rang);
    
    if(initiateur){
        printf("p%d : je suis initiateur.\n",id);
    }else{
        printf("p%d : je ne suis pas initiateur.\n",id);
    }
}



























/*
Algorithme de Hirschberg & Sinclair 
*/

/*
Fonction pow de math ne fonctionnait pas donc on a fait nôtre propre fonction power
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
Fonction qui envoye un message
*/

void SendMessage( int initiateur, int distance, int recepteur, int tag){
    msg[INIT] = initiateur;
    msg[DIST] = distance;
    //printf("p%d : j'envoie <init,dist> = <%d,%d> à %d avec le tag %d.\n",id,msg[INIT], msg[DIST], recepteur, tag);
    MPI_Send(msg, 2, MPI_INT, recepteur, tag, MPI_COMM_WORLD);
}

/*
Fonction qui initie l'étape 
*/

void InitierEtape(){
    cpt_in = 0;
    cpt_out = 0;
    int distance = power(2, etape);
    //printf("p%d : j'initie l'étape %d.\n",id, etape);
    SendMessage(id, distance, voisins[PREV], TAG_OUT);
    SendMessage(id, distance, voisins[NEXT], TAG_OUT);
    ++etape;
}

/*
Fonction qui reçoit un message de tag in
*/

void ReceiveIn(int emetteur, int initiateur, int distance){
    if( initiateur == id ){
        ++cpt_in;
        if(cpt_in == 2)     
            InitierEtape();
    }else{
        if(emetteur == voisins[PREV]){
            SendMessage(initiateur, -1, voisins[NEXT], TAG_IN);
        }else{
            SendMessage(initiateur, -1, voisins[PREV], TAG_IN);
        }

    }
}

/*
Fonction qui reçoit un message de tag out
*/

void ReceiveOut(int emetteur, int initiateur, int distance){
    if(initiateur_election == 0 || initiateur > id){
        state = BATTU;
        if(distance > 1){
            if(emetteur == voisins[NEXT]){
                SendMessage(initiateur, distance - 1, voisins[PREV], TAG_OUT);
            }else{
                SendMessage(initiateur, distance - 1, voisins[NEXT], TAG_OUT);
            }
        }else{
            if(emetteur == voisins[NEXT]){
                SendMessage(initiateur, -1, voisins[PREV], TAG_IN);
            }else{
                SendMessage(initiateur, -1, voisins[NEXT], TAG_IN);
            }
        }
    }else{
        if(initiateur == id && cpt_out == 0){
            ++cpt_out;
            state = ELU;
            taille = power(2, etape-1) - msg[DIST] + 1;
            printf("p%d : JE SUIS LE LEADER DE L'ANNEAU !!!\n",id);
            printf("p%d : la taille de l'anneau est %d.\n",id, taille);
            //SendMessage( id, 0, voisins[NEXT], TAG_LEADER);
            SendMessage( id, taille, voisins[NEXT], TAG_LEADER);
        }
    }
}

void ElectionLeader(MPI_Status status){

    // Initialisation des valeurs de départ
    initiateur_election = (initiateur == 1);
    
    if(initiateur_election){
        InitierEtape();
    }

    do{
        MPI_Recv(msg, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        //printf("p%d : reçu <init,dist> = <%d,%d> de %d.\n", id, msg[INIT], msg[DIST], status.MPI_SOURCE);

        switch (status.MPI_TAG) {
            case TAG_IN:
                ReceiveIn(status.MPI_SOURCE, msg[INIT], msg[DIST]);
                break;

            case TAG_OUT:
                ReceiveOut(status.MPI_SOURCE, msg[INIT], msg[DIST]);
                break;

            case TAG_LEADER:
                if(msg[INIT] == id){
                    // on ne fait rien
                }else{
                    printf("p%d : un leader a été élu -> p%d.\n",id,msg[INIT]);
                    SendMessage( msg[INIT], msg[DIST], voisins[NEXT], TAG_LEADER);
                    taille = msg[DIST];
                    printf("p%d : la taille de l'anneau est de %d.\n",id, taille);
                }
                break;

            default:
                printf("ERREUR.\n");

        }

    }while(status.MPI_TAG != TAG_LEADER);
}























/* Fonction transmission du tableaux d'id CHORD */

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
Le leader de l'anneau initie un tableau d'id CHORD et le transmet à son sucesseur.
Les autres noeuds de l'anneau transmettent juste le tableau.
Au second tour, tous les noeuds auront reçu l'ensemble des id chords de l'anneau.
*/

void TransmissionIDChord(int rang, MPI_Status status){
    tab_id = (int *) malloc( (taille+1) * sizeof(int));

    // Dans le cas où l'on est le leader de l'anneau
    if(state == ELU){
        /*
        tab_id_chord est le tableau des id chord.
        Chaque noeud va ajouter son id chord à l'emplacement correspondant à son rang MPI tq:
            tab_id_chord[0] = -1 pour que quant on l'ordonne avec un quicksort la valeur -1 reste à la position 0
            tab_id_chord[RANG_MPI] = id_chord
        
        */
        int tab_id_chord[taille + 1];

        tab_id_chord[0] = -1;

        for(int i = 1; i < taille + 1; ++i){
            tab_id_chord[i] = 0;
        }

        tab_id_chord[rang] = id;

        // Envoie le tableau des id chord à compléter
        MPI_Send(tab_id_chord, taille + 1, MPI_INT, voisins[NEXT], TAG_TAB_ID, MPI_COMM_WORLD);

        // Reçoit le tableau des id chord compléter
        MPI_Recv(tab_id, taille + 1, MPI_INT, voisins[PREV], TAG_TAB_ID, MPI_COMM_WORLD, &status);

        // Ordonne le tableau avant de l'envoyer
        qsort(tab_id, taille + 1, sizeof(int),fonction_cmp);

        // Envoie le tableau des id chord complet
        MPI_Send(tab_id, taille + 1, MPI_INT, voisins[NEXT], TAG_END, MPI_COMM_WORLD);

    }else{ 
        // Dans le cas où l'on n'est pas le leader de l'anneau

        // Reçoit le tableau des id chord à compléter
        MPI_Recv(tab_id, taille + 1, MPI_INT, voisins[PREV], TAG_TAB_ID, MPI_COMM_WORLD, &status);

        // Ajoute son id chord
        tab_id[rang] = id;

        // Envoie le tableau des id chord à compléter
        MPI_Send(tab_id, taille + 1, MPI_INT, voisins[NEXT], TAG_TAB_ID, MPI_COMM_WORLD);

        // Reçoit le tableau des id chord complet
        MPI_Recv(tab_id, taille + 1, MPI_INT, voisins[PREV], TAG_END, MPI_COMM_WORLD, &status);

        // Envoie le tableau des id chord complet
        MPI_Send(tab_id, taille + 1, MPI_INT, voisins[NEXT], TAG_END, MPI_COMM_WORLD);
    }
}

void CalculFingerTable(){
    int max_valeur = power(2, M);
    
    // Pour chaque finger
    for (int i = 0; i < M; ++i) {
        int valeur = (id + power(2, i)) % max_valeur;     // (id_chord + 2^i) mod 2^M
        int cpt_id = 1;                                      // Compteur du tableau des id CHORD

        // On cherche le plus petit id_chord plus grand que valeur
        while(cpt_id < NB_SITE && valeur > tab_id[cpt_id]){
            ++cpt_id;
        }

        if(cpt_id == NB_SITE) cpt_id = 1;                   // Si cpt_id == NB_SITE, alors valeur est compris entre le plus grand id_chord est 0, donc son responsable est le plus petit id_chord

        // On a trouvé le plus petit id_chord plus grand que valeur
        finger[i] = tab_id[cpt_id];
    }

    printf("p%d : Finger table avec M = %d.\n",id, M);
    printf("p%d : i\t id chord supérieur\n", id);
    for(int cpt = 0; cpt < M; ++cpt){
        printf("p%d : %d\t %d.\n",id, cpt, finger[cpt]);
    }
    printf("\n");
    free(tab_id);
}

























/* Fonction main */

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
        simulateur(rang);
    } else {
        Instanciation(rang, status);
        printf("\n");
        ElectionLeader(status);
        printf("\n");
        TransmissionIDChord(rang, status);
        printf("\n");
        CalculFingerTable();
    }
  
    MPI_Finalize();
    return 0;
}
