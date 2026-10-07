/* verify.c — exhaustive validation of ripes.c's IDA* solver against a
 * ground-truth BFS table. Host-only. ripes.c is included rather than copied,
 * so what is checked is the solver itself. Build:
 *
 *   cc -O2 -o verify verify.c
 *   ./verify
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Nodes ripes.c's search() has visited, counted through its COUNT_NODE hook. */
static uint64_t nodes;
#define COUNT_NODE() (++nodes)

#define main ripes_main
#include "ripes.c"
#undef main

enum {
    PERMUTATIONS = 5040,
    ORIENTATIONS = 729,
    STATES = PERMUTATIONS * ORIENTATIONS
};

/* ---- ground truth: BFS over the full product space ---------------------- */

static uint8_t *build_truth(uint8_t *diameter)
{
    uint8_t *dist = malloc(STATES);
    uint32_t *queue = malloc((size_t) STATES * sizeof(uint32_t));
    uint32_t head = 0, tail = 0;
    uint8_t max = 0;

    if (!dist || !queue) {
        free(dist);
        free(queue);
        return NULL;
    }
    memset(dist, 0xFF, STATES);

    dist[0] = 0;              /* solved = permutation rank 0, orientation 0 */
    queue[tail++] = 0;

    while (head < tail) {
        uint32_t rank = queue[head++];
        uint16_t p = (uint16_t) (rank / ORIENTATIONS);
        uint16_t o = (uint16_t) (rank % ORIENTATIONS);
        uint8_t d = dist[rank];

        for (uint8_t face = 0; face < FACES; ++face) {
            uint16_t next_p = p, next_o = o;
            for (uint8_t turn = 0; turn < 3; ++turn) {
                uint32_t next;
                next_p = permutation_table[next_p][face] / ROW_BYTES;
                next_o = orientation_table[next_o][face] / ROW_BYTES;
                next = (uint32_t) next_p * ORIENTATIONS + next_o;
                if (dist[next] != 0xFF)
                    continue;
                dist[next] = (uint8_t) (d + 1);
                if (d + 1 > max)
                    max = (uint8_t) (d + 1);
                queue[tail++] = next;
            }
        }
    }
    free(queue);
    *diameter = max;
    return (tail == STATES) ? dist : (free(dist), NULL);
}

/* ---- enumerating states -----------------------------------------------
 *
 * Turns a rank of the BFS table back into a state for solve(). This must agree
 * with solver.c's rank_state, which ripes_tables.h is indexed by.
 */
static void unrank(uint32_t rank, state_t *s)
{
    uint8_t available[CUBIES];
    uint32_t pr = rank / ORIENTATIONS, or = rank % ORIENTATIONS;
    uint32_t f = 720;
    uint8_t sum = 0;

    for (uint8_t i = 0; i < CUBIES; ++i)
        available[i] = i;
    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint32_t q = f ? pr / f : 0;
        s->p[i] = available[q];
        for (uint8_t j = (uint8_t) q; j + 1 < CUBIES - i; ++j)
            available[j] = available[j + 1];
        if (f) {
            pr %= f;
            f /= (CUBIES - 1 - i) ? (CUBIES - 1 - i) : 1;
        }
    }
    for (int8_t i = CUBIES - 2; i >= 0; --i) {
        s->o[i] = (uint8_t) (or % 3);
        or /= 3;
        sum = (uint8_t) (sum + s->o[i]);
    }
    s->o[CUBIES - 1] = (uint8_t) ((3U - sum % 3U) % 3U);
}

/* ---- checks ------------------------------------------------------------- */

int main(void)
{
    uint8_t diameter, *truth;
    clock_t t0;

    printf("building ground-truth BFS table...\n");
    t0 = clock();
    truth = build_truth(&diameter);
    if (!truth) {
        fprintf(stderr, "ground truth build failed\n");
        return 1;
    }
    printf("  %d states, diameter %u, %.3f s\n\n", STATES, diameter,
           (double) (clock() - t0) / CLOCKS_PER_SEC);

    /* ---- H2: the pruning tables are complete and correctly bounded ---- */
    /* Each successor must also be the byte offset of a row, since ripes.c
     * adds column offsets to it and reads there unchecked. */
    {
        uint16_t pmax = 0, omax = 0;
        int ok = 1;
        for (uint16_t i = 0; i < PERMUTATIONS; ++i) {
            if (permutation_table[i][FACES] == 0xFF) {
                printf("H2 FAIL: permutation depth %u unfilled\n", i);
                ok = 0;
            }
            if (permutation_table[i][FACES] > pmax)
                pmax = permutation_table[i][FACES];
            for (uint8_t face = 0; face < FACES; ++face)
                if (permutation_table[i][face] % ROW_BYTES ||
                    permutation_table[i][face] >= PERMUTATIONS * ROW_BYTES) {
                    printf("H2 FAIL: permutation_table[%u][%u] = %u\n", i,
                           face, permutation_table[i][face]);
                    ok = 0;
                }
        }
        for (uint16_t i = 0; i < ORIENTATIONS; ++i) {
            if (orientation_table[i][FACES] == 0xFF) {
                printf("H2 FAIL: orientation depth %u unfilled\n", i);
                ok = 0;
            }
            if (orientation_table[i][FACES] > omax)
                omax = orientation_table[i][FACES];
            for (uint8_t face = 0; face < FACES; ++face)
                if (orientation_table[i][face] % ROW_BYTES ||
                    orientation_table[i][face] >= ORIENTATIONS * ROW_BYTES) {
                    printf("H2 FAIL: orientation_table[%u][%u] = %u\n", i,
                           face, orientation_table[i][face]);
                    ok = 0;
                }
        }
        if (permutation_table[0][FACES] != 0) {
            printf("H2 FAIL: permutation depth of solved = %u\n",
                   permutation_table[0][FACES]);
            ok = 0;
        }
        if (orientation_table[0][FACES] != 0) {
            printf("H2 FAIL: orientation depth of solved = %u\n",
                   orientation_table[0][FACES]);
            ok = 0;
        }
        printf("H2 %s: tables filled, max depths %u / %u, solved entries 0\n",
               ok ? "PASS" : "FAIL", pmax, omax);
    }

    /* ---- H1: both heuristics are admissible on every state ---- */
    t0 = clock();
    {
        uint64_t violations = 0;
        uint8_t worst_slack = 0;
        for (uint32_t rank = 0; rank < STATES; ++rank) {
            uint16_t p = (uint16_t) (rank / ORIENTATIONS);
            uint16_t o = (uint16_t) (rank % ORIENTATIONS);
            uint8_t d = truth[rank];
            uint8_t h = (uint8_t) (permutation_table[p][FACES] >
                                           orientation_table[o][FACES]
                                       ? permutation_table[p][FACES]
                                       : orientation_table[o][FACES]);
            if (h > d) {
                if (violations < 5)
                    printf("H1 FAIL: rank %u, h=%u > d=%u\n", rank, h, d);
                ++violations;
            } else if ((uint8_t) (d - h) > worst_slack) {
                worst_slack = (uint8_t) (d - h);
            }
        }
        printf("H1 %s: %llu violations over %d states, worst slack %u (%.3f s)\n",
               violations ? "FAIL" : "PASS", (unsigned long long) violations,
               STATES, worst_slack, (double) (clock() - t0) / CLOCKS_PER_SEC);
    }

    /* ---- H3: the solver returns the exact distance on every state ---- */
    /* ---- plus node counts, and the distance-11 worst cases ---- */
    t0 = clock();
    {
        uint64_t mismatches = 0, total_nodes = 0, d11_nodes = 0;
        uint64_t max_nodes = 0;
        uint32_t max_rank = 0, d11_count = 0;
        /* top five distance-11 states by node count */
        uint64_t top_n[5] = {0};
        uint32_t top_r[5] = {0};

        for (uint32_t rank = 0; rank < STATES; ++rank) {
            state_t s;
            int length;

            unrank(rank, &s);
            nodes = 0;
            length = solve(&s);
            total_nodes += nodes;

            if (length != truth[rank]) {
                if (mismatches < 5)
                    printf("H3 FAIL: rank %u, got %d, expected %u\n", rank,
                           length, truth[rank]);
                ++mismatches;
            }
            if (nodes > max_nodes) {
                max_nodes = nodes;
                max_rank = rank;
            }
            if (truth[rank] == 11) {
                ++d11_count;
                d11_nodes += nodes;
                for (int k = 0; k < 5; ++k) {
                    if (nodes > top_n[k]) {
                        for (int j = 4; j > k; --j) {
                            top_n[j] = top_n[j - 1];
                            top_r[j] = top_r[j - 1];
                        }
                        top_n[k] = nodes;
                        top_r[k] = rank;
                        break;
                    }
                }
            }
            if ((rank & 0xFFFFF) == 0)
                printf("  ... %u / %d\n", rank, STATES);
        }

        printf("\nH3 %s: %llu mismatches over %d states (%.1f s)\n",
               mismatches ? "FAIL" : "PASS", (unsigned long long) mismatches,
               STATES, (double) (clock() - t0) / CLOCKS_PER_SEC);
        printf("  mean nodes      %.0f\n", (double) total_nodes / STATES);
        printf("  max nodes       %llu at rank %u (distance %u)\n",
               (unsigned long long) max_nodes, max_rank, truth[max_rank]);
        printf("  distance-11:    %u states, mean %.0f nodes\n", d11_count,
               (double) d11_nodes / d11_count);
        printf("  worst five distance-11 states:\n");
        for (int k = 0; k < 5; ++k) {
            state_t s;
            unrank(top_r[k], &s);
            printf("    ");
            for (int i = 0; i < CUBIES; ++i)
                printf("%u", s.p[i] + 1);
            for (int i = 0; i < CUBIES; ++i)
                printf("%u", s.o[i] + 1);
            printf("  rank %-8u %llu nodes\n", top_r[k],
                   (unsigned long long) top_n[k]);
        }
    }

    /* ---- H4 ---- */
    printf("\nH4 N/A: no packed accessor — the tables are plain arrays.\n");

    free(truth);
    return 0;
}