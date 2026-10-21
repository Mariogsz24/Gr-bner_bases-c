#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>


#define MAX_VARS 2
#define MAX_TERMS 1000
#define MAX_POLIS 1000
#define EPSILON 1e-6

int selected_order = 1;

typedef struct {
    double coeff;
    int exps[MAX_VARS];
} Term;

typedef struct {
    Term terms[MAX_TERMS];
    int num_terms;
} Poli;

typedef struct {
    Poli polis[MAX_POLIS];
    int count;
} PoliList;

int cmp_terms(const void *a, const void *b) {
    Term *t1 = (Term *)a;
    Term *t2 = (Term *)b;

    if (selected_order == 1) {
        for (int i = 0; i < MAX_VARS; i++) {
            if (t1->exps[i] != t2->exps[i]) {
                return t2->exps[i] - t1->exps[i];
            }
        }
        return 0;
    } 
    else if (selected_order == 2) {
        int deg1 = 0, deg2 = 0;
        for (int i = 0; i < MAX_VARS; i++) {
            deg1 += t1->exps[i];
            deg2 += t2->exps[i];
        }
        if (deg1 != deg2) {
            return deg2 - deg1;
       }
        
        for (int i = 0; i < MAX_VARS; i++) {
            if (t1->exps[i] != t2->exps[i]) {
                return t2->exps[i] - t1->exps[i];
            }
        }
        return 0;
    } 
    else if (selected_order == 3) {
        int deg1 = 0, deg2 = 0;
        for (int i = 0; i < MAX_VARS; i++) {
            deg1 += t1->exps[i];
            deg2 += t2->exps[i];
        }
        if (deg1 != deg2) {
            return deg2 - deg1; 
        }
        
        for (int i = MAX_VARS - 1; i >= 0; i--) {
            if (t1->exps[i] != t2->exps[i]) {
                return t1->exps[i] - t2->exps[i]; 
            }
        }
        return 0;
    }
    return 0;
}

void simplify_poli(Poli *p) {
    if (p->num_terms == 0) return;
    
    qsort(p->terms, p->num_terms, sizeof(Term), cmp_terms);
    
    Poli temp;
    temp.num_terms = 0;
    
    for (int i = 0; i < p->num_terms; i++) {
        if (fabs(p->terms[i].coeff) < EPSILON) continue;
        
        if (temp.num_terms > 0 && cmp_terms(&p->terms[i], &temp.terms[temp.num_terms - 1]) == 0) {
            temp.terms[temp.num_terms - 1].coeff += p->terms[i].coeff;
            if (fabs(temp.terms[temp.num_terms - 1].coeff) < EPSILON) {
                temp.num_terms--; 
            }
        } else {
            temp.terms[temp.num_terms++] = p->terms[i];
        }
    }
    *p = temp;
}

Poli sumar_poli(Poli p1, Poli p2) {
    Poli res; res.num_terms = 0;
    for (int i = 0; i < p1.num_terms; i++) res.terms[res.num_terms++] = p1.terms[i];
    for (int i = 0; i < p2.num_terms; i++) res.terms[res.num_terms++] = p2.terms[i];
    simplify_poli(&res);
    return res;
}

Poli mult_term_poli(Poli p, Term t) {
    Poli res; res.num_terms = p.num_terms;
    for (int i = 0; i < p.num_terms; i++) {
        res.terms[i].coeff = p.terms[i].coeff * t.coeff;
        for (int v = 0; v < MAX_VARS; v++) res.terms[i].exps[v] = p.terms[i].exps[v] + t.exps[v];
    }
    simplify_poli(&res);
    return res;
}

Poli restar_poli(Poli p1, Poli p2) {
    Term negativo = { -1.0, {0} };
    return sumar_poli(p1, mult_term_poli(p2, negativo));
}

bool divides(Term t1, Term t2) {
    for (int i = 0; i < MAX_VARS; i++) if (t1.exps[i] < t2.exps[i]) return false;
    return true;
}

Term div_term(Term t1, Term t2) {
    Term res; res.coeff = t1.coeff / t2.coeff;
    for (int i = 0; i < MAX_VARS; i++) res.exps[i] = t1.exps[i] - t2.exps[i];
    return res;
}

bool son_coprimos(Term t1, Term t2) {
    for (int i = 0; i < MAX_VARS; i++) {
        
        if (t1.exps[i] > 0 && t2.exps[i] > 0) {
            return false; 
        }
    }
    return true; 
}

Term lcm_term(Term t1, Term t2) {
    Term res = {1.0, {0}};
    for (int i = 0; i < MAX_VARS; i++) res.exps[i] = (t1.exps[i] > t2.exps[i]) ? t1.exps[i] : t2.exps[i];
    return res;
}

Poli division_multivariante(Poli p, PoliList *F) {
    Poli r; r.num_terms = 0;
    Poli temp_p = p;

    while (temp_p.num_terms > 0) {
        bool dividido = false;
        for (int i = 0; i < F->count; i++) {
            if (F->polis[i].num_terms > 0 && divides(temp_p.terms[0], F->polis[i].terms[0])) {
                Term qt = div_term(temp_p.terms[0], F->polis[i].terms[0]);
                Poli subtraido = mult_term_poli(F->polis[i], qt);
                temp_p = restar_poli(temp_p, subtraido);
                dividido = true;
                break;
            }
        }
        if (!dividido) {
            r.terms[r.num_terms++] = temp_p.terms[0];
            Poli next_p; next_p.num_terms = 0;
            for(int k = 1; k < temp_p.num_terms; k++) next_p.terms[next_p.num_terms++] = temp_p.terms[k];
            temp_p = next_p;
        }
    }
    simplify_poli(&r);
    return r;
}

Poli s_poli(Poli p1, Poli p2) {
    Term lcm = lcm_term(p1.terms[0], p2.terms[0]);
    Term t1 = div_term(lcm, p1.terms[0]);
    Term t2 = div_term(lcm, p2.terms[0]);
    return restar_poli(mult_term_poli(p1, t1), mult_term_poli(p2, t2));
}

void buchberger(PoliList *G) {
    bool added = true;
    while (added) {
        added = false;
        int current_count = G->count;
        for (int i = 0; i < current_count; i++) {
            for (int j = i + 1; j < current_count; j++) {

		if (son_coprimos(G->polis[i].terms[0], G->polis[j].terms[0])) {
                    continue; 
                }

                Poli sp = s_poli(G->polis[i], G->polis[j]);
                if (sp.num_terms == 0) continue;
                
                Poli rem = division_multivariante(sp, G);
                if (rem.num_terms > 0) {
                    G->polis[G->count++] = rem;
                    added = true;
                }
            }
        }
    }
}

void make_monic(Poli *p) {
    if (p->num_terms == 0) return;
    double lead_coeff = p->terms[0].coeff;
    for (int i = 0; i < p->num_terms; i++) p->terms[i].coeff /= lead_coeff;
}

void grobner_reducido(PoliList *G) {
    buchberger(G); 

    bool active[MAX_POLIS];
    for (int i = 0; i < G->count; i++) active[i] = true;

    for (int i = 0; i < G->count; i++) {
        for (int j = 0; j < G->count; j++) {
            if (i != j && active[i] && active[j]) {
                if (divides(G->polis[i].terms[0], G->polis[j].terms[0])) {
                    active[i] = false; 
                    break;
                }
            }
        }
    }

    PoliList *minimalG = (PoliList *)malloc(sizeof(PoliList));
    if (minimalG == NULL) exit(1);
    minimalG->count = 0;

    for (int i = 0; i < G->count; i++) {
        if (active[i]) {
            make_monic(&G->polis[i]);
            minimalG->polis[minimalG->count++] = G->polis[i];
        }
    }

    *G = *minimalG;
    free(minimalG);

    for (int i = 0; i < G->count; i++) {
        Poli p = G->polis[i];
        
        PoliList *G_except_p = (PoliList *)malloc(sizeof(PoliList));
        if (G_except_p == NULL) exit(1);
        G_except_p->count = 0;
        
        for (int j = 0; j < G->count; j++) {
            if (i != j) G_except_p->polis[G_except_p->count++] = G->polis[j];
        }
        
        G->polis[i] = division_multivariante(p, G_except_p);
        make_monic(&G->polis[i]); 
        free(G_except_p);
    }
}

void print_poli(Poli p) {
    if (p.num_terms == 0) { printf("0\n"); return; }
    const char var_names[] = "xyztuvwqrsabcdefghijklmnop";
    
    for (int i = 0; i < p.num_terms; i++) {
        if (i > 0 && p.terms[i].coeff > 0) printf(" + ");
        else if (p.terms[i].coeff < 0) printf(" - ");
        
        double c = fabs(p.terms[i].coeff);

        int ceros = 1;

        for (int j = 0; j < MAX_VARS; j++) {
            if (p.terms[i].exps[j] != 0) {
                ceros = 0;
                break;
            }
        }

        if (fabs(c - 1.0) > EPSILON || ceros ) {
            printf("%.2f", c);
        }

        for (int j = 0; j < MAX_VARS; j++) {
            if (p.terms[i].exps[j] > 0) {
                
                char var_char = (j < sizeof(var_names) - 1) ? var_names[j] : '?'; 
                
                printf("%c^%d", var_char, p.terms[i].exps[j]);
            }
        }
    }
    printf("\n");
}

int main() {
    
    printf("   BASE DE GRÖBNER REDUCIDA \n");
    printf("Orden monomial deseado:\n");
    printf("  1 - LEX \n");
    printf("  2 - GRLEX \n");
    printf("  3 - GREVLEX \n");
    printf("Ingrese su opcion: ");
    
    if (scanf("%d", &selected_order) != 1 || selected_order < 1 || selected_order > 3) {
        printf("Opcion invalida. Se utilizara LEX (1) por defecto.\n");
        selected_order = 1;
    }

    PoliList *G = (PoliList *)malloc(sizeof(PoliList));
    if (G == NULL) exit(1);
    G->count = 0;
    
    // Definicion de los polinomios (coeficiente, {Exponentes de las variables})
    Poli f1; f1.num_terms = 3;    
    f1.terms[0] = (Term){3.0, {5, 0}};
    f1.terms[1] = (Term){4.0, {0, 5}};
    f1.terms[2] = (Term){-7.0, {0, 0}};

    
    Poli f2; f2.num_terms = 3;    
    f2.terms[0] = (Term){2.0, {2, 1}};
    f2.terms[1] = (Term){3.0, {0, 3}};
    f2.terms[2] = (Term){-1.0, {0, 0}};

   
    simplify_poli(&f1); //Ordenacion monomial
    simplify_poli(&f2);
    

    G->polis[G->count++] = f1; //Creacion de la base inicial
    G->polis[G->count++] = f2;
    

    printf("\n Polinomios Iniciales: \n");
    for (int i = 0; i < G->count; i++) {
        printf(" f%d = ", i + 1);
        print_poli(G->polis[i]);
    }

    clock_t start_time = clock();
    grobner_reducido(G);
    clock_t end_time = clock();

    char *names[] = {"", "LEX", "GRLEX", "GREVLEX"};
    printf("\n Base de Gröbner Minimal Reducida con orden %s \n", names[selected_order]);
    for (int i = 0; i < G->count; i++) {
        printf(" g%d = ", i + 1);
        print_poli(G->polis[i]);
    }
    
    double time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Tiempo de ejecución : %.6f segundos.\n", time);
    
    free(G);
    
    return 0;
}