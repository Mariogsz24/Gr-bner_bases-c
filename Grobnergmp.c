#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <gmp.h> //librería GMP

#define MAX_VARS 2
#define MAX_TERMS 1000
#define MAX_POLIS 1000

int selected_order = 1;

typedef struct {
    mpq_t coeff;
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

//inicializar un polinomio (necesario en GMP para reservar memoria)
void init_poli(Poli *p) {
    p->num_terms = 0;
}

//copiar un término
void copy_term(Term *dest, const Term *src) {
    mpq_init(dest->coeff);
    mpq_set(dest->coeff, src->coeff);
    for (int i = 0; i < MAX_VARS; i++) dest->exps[i] = src->exps[i];
}

//copiar un polinomio
void copy_poli(Poli *dest, const Poli *src) {
    dest->num_terms = src->num_terms;
    for (int i = 0; i < src->num_terms; i++) {
        copy_term(&dest->terms[i], &src->terms[i]);
    }
}

// Liberar memoria de los coeficientes de un polinomio
void clear_poli(Poli *p) {
    for (int i = 0; i < p->num_terms; i++) {
        mpq_clear(p->terms[i].coeff);
    }
    p->num_terms = 0;
}

int cmp_terms(const void *a, const void *b) {
    Term *t1 = (Term *)a;
    Term *t2 = (Term *)b;

    if (selected_order == 1) {
        for (int i = 0; i < MAX_VARS; i++) {
            if (t1->exps[i] != t2->exps[i]) return t2->exps[i] - t1->exps[i];
        }
    } else if (selected_order == 2) {
        int deg1 = 0, deg2 = 0;
        for (int i = 0; i < MAX_VARS; i++) { deg1 += t1->exps[i]; deg2 += t2->exps[i]; }
        if (deg1 != deg2) return deg2 - deg1;
        for (int i = 0; i < MAX_VARS; i++) {
            if (t1->exps[i] != t2->exps[i]) return t2->exps[i] - t1->exps[i];
        }
    } else if (selected_order == 3) {
        int deg1 = 0, deg2 = 0;
        for (int i = 0; i < MAX_VARS; i++) { deg1 += t1->exps[i]; deg2 += t2->exps[i]; }
        if (deg1 != deg2) return deg2 - deg1;
        for (int i = MAX_VARS - 1; i >= 0; i--) {
            if (t1->exps[i] != t2->exps[i]) return t1->exps[i] - t2->exps[i];
        }
    }
    return 0;
}

void simplify_poli(Poli *p) {
    if (p->num_terms == 0) return;
    
    
    qsort(p->terms, p->num_terms, sizeof(Term), cmp_terms);
    
    Poli temp;
    init_poli(&temp);
    
    for (int i = 0; i < p->num_terms; i++) {
        if (mpq_sgn(p->terms[i].coeff) == 0) continue; 
        
        if (temp.num_terms > 0 && cmp_terms(&p->terms[i], &temp.terms[temp.num_terms - 1]) == 0) {
            mpq_add(temp.terms[temp.num_terms - 1].coeff, temp.terms[temp.num_terms - 1].coeff, p->terms[i].coeff);
            if (mpq_sgn(temp.terms[temp.num_terms - 1].coeff) == 0) {
                mpq_clear(temp.terms[temp.num_terms - 1].coeff);
                temp.num_terms--; 
            }
        } else {
            copy_term(&temp.terms[temp.num_terms++], &p->terms[i]);
        }
    }
    clear_poli(p);
    copy_poli(p, &temp);
    clear_poli(&temp);
}

Poli sumar_poli(Poli p1, Poli p2) {
    Poli res; init_poli(&res);
    for (int i = 0; i < p1.num_terms; i++) copy_term(&res.terms[res.num_terms++], &p1.terms[i]);
    for (int i = 0; i < p2.num_terms; i++) copy_term(&res.terms[res.num_terms++], &p2.terms[i]);
    simplify_poli(&res);
    return res;
}

Poli mult_term_poli(Poli p, Term t) {
    Poli res; init_poli(&res);
    for (int i = 0; i < p.num_terms; i++) {
        mpq_init(res.terms[i].coeff);
        mpq_mul(res.terms[i].coeff, p.terms[i].coeff, t.coeff);
        for (int v = 0; v < MAX_VARS; v++) res.terms[i].exps[v] = p.terms[i].exps[v] + t.exps[v];
        res.num_terms++;
    }
    simplify_poli(&res);
    return res;
}

Poli restar_poli(Poli p1, Poli p2) {
    Term negativo;
    mpq_init(negativo.coeff);
    mpq_set_si(negativo.coeff, -1, 1);
    for (int i = 0; i < MAX_VARS; i++) negativo.exps[i] = 0;
    
    Poli temp = mult_term_poli(p2, negativo);
    Poli res = sumar_poli(p1, temp);
    
    clear_poli(&temp);
    mpq_clear(negativo.coeff);
    return res;
}

bool divides(Term t1, Term t2) {
    for (int i = 0; i < MAX_VARS; i++) if (t1.exps[i] < t2.exps[i]) return false;
    return true;
}

Term div_term(Term t1, Term t2) {
    Term res; 
    mpq_init(res.coeff);
    mpq_div(res.coeff, t1.coeff, t2.coeff);
    for (int i = 0; i < MAX_VARS; i++) res.exps[i] = t1.exps[i] - t2.exps[i];
    return res;
}

bool son_coprimos(Term t1, Term t2) {
    for (int i = 0; i < MAX_VARS; i++) {
        if (t1.exps[i] > 0 && t2.exps[i] > 0) return false; 
    }
    return true; 
}

Term lcm_term(Term t1, Term t2) {
    Term res;
    mpq_init(res.coeff);
    mpq_set_si(res.coeff, 1, 1); // LCM coefficient es 1
    for (int i = 0; i < MAX_VARS; i++) res.exps[i] = (t1.exps[i] > t2.exps[i]) ? t1.exps[i] : t2.exps[i];
    return res;
}

Poli division_multivariante(Poli p, PoliList *F) {
    Poli r; init_poli(&r);
    Poli temp_p; init_poli(&temp_p);
    copy_poli(&temp_p, &p);

    while (temp_p.num_terms > 0) {
        bool dividido = false;
        for (int i = 0; i < F->count; i++) {
            if (F->polis[i].num_terms > 0 && divides(temp_p.terms[0], F->polis[i].terms[0])) {
                Term qt = div_term(temp_p.terms[0], F->polis[i].terms[0]);
                Poli subtraido = mult_term_poli(F->polis[i], qt);
                
                Poli nueva_temp = restar_poli(temp_p, subtraido);
                clear_poli(&temp_p);
                temp_p = nueva_temp;
                
                clear_poli(&subtraido);
                mpq_clear(qt.coeff);
                dividido = true;
                break;
            }
        }
        if (!dividido) {
            copy_term(&r.terms[r.num_terms++], &temp_p.terms[0]);
            Poli next_p; init_poli(&next_p);
            for(int k = 1; k < temp_p.num_terms; k++) copy_term(&next_p.terms[next_p.num_terms++], &temp_p.terms[k]);
            clear_poli(&temp_p);
            temp_p = next_p;
        }
    }
    clear_poli(&temp_p);
    simplify_poli(&r);
    return r;
}

Poli s_poli(Poli p1, Poli p2) {
    Term lcm = lcm_term(p1.terms[0], p2.terms[0]);
    Term t1 = div_term(lcm, p1.terms[0]);
    Term t2 = div_term(lcm, p2.terms[0]);
    
    Poli m1 = mult_term_poli(p1, t1);
    Poli m2 = mult_term_poli(p2, t2);
    Poli res = restar_poli(m1, m2);
    
    mpq_clear(lcm.coeff);
    mpq_clear(t1.coeff);
    mpq_clear(t2.coeff);
    clear_poli(&m1);
    clear_poli(&m2);
    
    return res;
}

void buchberger(PoliList *G) {
    bool added = true;
    while (added) {
        added = false;
        int current_count = G->count;
        for (int i = 0; i < current_count; i++) {
            for (int j = i + 1; j < current_count; j++) {
		        if (son_coprimos(G->polis[i].terms[0], G->polis[j].terms[0])) continue; 

                Poli sp = s_poli(G->polis[i], G->polis[j]);
                if (sp.num_terms == 0) {
                    clear_poli(&sp);
                    continue;
                }
                
                Poli rem = division_multivariante(sp, G);
                if (rem.num_terms > 0) {
                    init_poli(&G->polis[G->count]);
                    copy_poli(&G->polis[G->count++], &rem);
                    added = true;
                }
                clear_poli(&sp);
                clear_poli(&rem);
            }
        }
    }
}

void make_monic(Poli *p) {
    if (p->num_terms == 0) return;
    
    mpq_t lead_coeff;
    mpq_init(lead_coeff);
    mpq_set(lead_coeff, p->terms[0].coeff);
    
    for (int i = 0; i < p->num_terms; i++) {
        mpq_div(p->terms[i].coeff, p->terms[i].coeff, lead_coeff);
    }
    mpq_clear(lead_coeff);
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
            init_poli(&minimalG->polis[minimalG->count]);
            copy_poli(&minimalG->polis[minimalG->count++], &G->polis[i]);
        }
    }

    for(int i = 0; i < G->count; i++) clear_poli(&G->polis[i]);
    
    G->count = minimalG->count;
    for (int i = 0; i < minimalG->count; i++) {
        init_poli(&G->polis[i]);
        copy_poli(&G->polis[i], &minimalG->polis[i]);
        clear_poli(&minimalG->polis[i]);
    }
    free(minimalG);

    for (int i = 0; i < G->count; i++) {
        Poli p; init_poli(&p);
        copy_poli(&p, &G->polis[i]);
        
        PoliList *G_except_p = (PoliList *)malloc(sizeof(PoliList));
        if (G_except_p == NULL) exit(1);
        G_except_p->count = 0;
        
        for (int j = 0; j < G->count; j++) {
            if (i != j) {
                init_poli(&G_except_p->polis[G_except_p->count]);
                copy_poli(&G_except_p->polis[G_except_p->count++], &G->polis[j]);
            }
        }
        
        clear_poli(&G->polis[i]);
        G->polis[i] = division_multivariante(p, G_except_p);
        make_monic(&G->polis[i]); 
        
        clear_poli(&p);
        for(int j = 0; j < G_except_p->count; j++) clear_poli(&G_except_p->polis[j]);
        free(G_except_p);
    }
}

void print_poli(Poli p) {
    if (p.num_terms == 0) { printf("0\n"); return; }
    const char var_names[] = "xyztuvwqrsabcdefghijklmnop";
    
    for (int i = 0; i < p.num_terms; i++) {
        int sgn = mpq_sgn(p.terms[i].coeff);
        
        if (i > 0 && sgn > 0) printf(" + ");
        else if (sgn < 0) {
            if(i == 0) printf("-");
            else printf(" - ");
        }
        
        mpq_t abs_coeff;
        mpq_init(abs_coeff);
        mpq_abs(abs_coeff, p.terms[i].coeff);

        int ceros = 1;
        for (int j = 0; j < MAX_VARS; j++) {
            if (p.terms[i].exps[j] != 0) {
                ceros = 0; break;
            }
        }

        // Si no es 1, o si es el término independiente, imprimimos el coeficiente
        if (mpq_cmp_si(abs_coeff, 1, 1) != 0 || ceros) {
            gmp_printf("%Qd", abs_coeff);
        }
        mpq_clear(abs_coeff);

        for (int j = 0; j < MAX_VARS; j++) {
            if (p.terms[i].exps[j] > 0) {
                char var_char = (j < sizeof(var_names) - 1) ? var_names[j] : '?'; 
                if(p.terms[i].exps[j] == 1)
                    printf("%c", var_char);
                else
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
    
    // Definicion manual de los polinomios con inicialización de GMP
    Poli f1; init_poli(&f1); f1.num_terms = 3;    
    mpq_init(f1.terms[0].coeff); mpq_set_si(f1.terms[0].coeff, 3, 1); f1.terms[0].exps[0] = 5; f1.terms[0].exps[1] = 0;
    mpq_init(f1.terms[1].coeff); mpq_set_si(f1.terms[1].coeff, 4, 1); f1.terms[1].exps[0] = 0; f1.terms[1].exps[1] = 5;
    mpq_init(f1.terms[2].coeff); mpq_set_si(f1.terms[2].coeff, -7, 1); f1.terms[2].exps[0] = 0; f1.terms[2].exps[1] = 0;
    
    Poli f2; init_poli(&f2); f2.num_terms = 3;    
    mpq_init(f2.terms[0].coeff); mpq_set_si(f2.terms[0].coeff, 2, 1); f2.terms[0].exps[0] = 2; f2.terms[0].exps[1] = 1;
    mpq_init(f2.terms[1].coeff); mpq_set_si(f2.terms[1].coeff, 3, 1); f2.terms[1].exps[0] = 0; f2.terms[1].exps[1] = 3;
    mpq_init(f2.terms[2].coeff); mpq_set_si(f2.terms[2].coeff, -1, 1); f2.terms[2].exps[0] = 0; f2.terms[2].exps[1] = 0;

    simplify_poli(&f1);
    simplify_poli(&f2);
    
    init_poli(&G->polis[0]); copy_poli(&G->polis[G->count++], &f1);
    init_poli(&G->polis[1]); copy_poli(&G->polis[G->count++], &f2);

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
    
    for (int i = 0; i < G->count; i++) clear_poli(&G->polis[i]);
    free(G);
    clear_poli(&f1);
    clear_poli(&f2);
    
    return 0;
}