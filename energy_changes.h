#ifndef ENERGY_CHANGES_H
#define ENERGY_CHANGES_H

#include <cmath>
#include "global.h"
using namespace Eigen;

bool accept_move(int thread_num, const double& delta_E) {
    return generators[thread_num].disReal() <= exp(-delta_E);
}

double delta_E_other(std::vector<Vector3i>& polym, std::vector<Vector3i>& other_polym, int site, Vector3i prop_move1, int thread_num, int lin_pol) {
    double energy_change = 0;
//  int stage = thread_num % number_of_stages;
//  int stage = thread_num; //GG: when length.size()=numofthreads (needed for fork distribution)
    int stage = thread_num / replicates_per_stage;

    //Energies from contacts if (site+1)%pol_length is a physical monomer
    if((site+1)%reduction_factor==0){
        int red_i_plus_one=((site + 1) % pol_length)/reduction_factor;
        if (locations[thread_num].find({ prop_move1[0],prop_move1[1],prop_move1[2] }) != locations[thread_num].end()) {
            for (auto elem : locations[thread_num][{ prop_move1[0],prop_move1[1],prop_move1[2] }]) {
                energy_change += Interaction_E[elem][red_i_plus_one];
            }
        }
        if (locations_lin[thread_num].find({ prop_move1[0],prop_move1[1],prop_move1[2] }) != locations_lin[thread_num].end()) {
            for (auto elem : locations_lin[thread_num][{ prop_move1[0],prop_move1[1],prop_move1[2] }]) {
                energy_change += Interaction_E[elem][red_i_plus_one];
            }
        }
        for (auto elem : locations[thread_num][{ polym[(site + 1) % pol_length][0],polym[(site + 1) % pol_length][1],polym[(site + 1) % pol_length][2] }]) {
            if (elem != red_i_plus_one) {
                energy_change -= Interaction_E[elem][red_i_plus_one];
            }
        }
        for (auto elem : locations_lin[thread_num][{ polym[(site + 1) % pol_length][0],polym[(site + 1) % pol_length][1],polym[(site + 1) % pol_length][2] }]) {
            if (elem != red_i_plus_one) {
                energy_change -= Interaction_E[elem][red_i_plus_one];
            }
        }
    }


    //Non-contact energy contributions
    if (lin_length[stage] !=0) {
        if (site + 1 == oriC && lin_pol) { // lin_polymer oriC
//            if (pol_close_pole(site+1,thread_num)) {
//                energy_change += beta[stage] * (std::abs(prop_move1[2]-pole[thread_num]) - std::abs(polym[site + 1][2] - pole[thread_num]));
//            }
//            else {
//                energy_change += alpha[stage] * (std::abs(prop_move1[2]-pole[thread_num]) - std::abs(polym[site + 1][2] - pole[thread_num]));
//            }
            energy_change += beta[stage] * (prop_move1[2]  - polym[site + 1][2] );
            energy_change += beta2[stage] *(pow(prop_move1[2] -xp_z_far_simunits[stage],2) - pow(polym[site + 1][2] -xp_z_far_simunits[stage],2));
        }
        else if (site + 1 == oriC && !lin_pol) {  // polymer oriC
//            if (!pol_close_pole(site+1,thread_num)) {
//                energy_change += beta[stage] * (std::abs(prop_move1[2]-pole[thread_num]) - std::abs(polym[site + 1][2] - pole[thread_num]));
//            }
//            else {
//                energy_change += alpha[stage] * (std::abs(prop_move1[2]-pole[thread_num]) - std::abs(polym[site + 1][2] - pole[thread_num]));
//            }
            energy_change += alpha[stage] * (prop_move1[2] - polym[site + 1][2] );
            energy_change += alpha2[stage] * (pow(prop_move1[2] -xp_z_close_simunits[stage],2) - pow(polym[site + 1][2] -xp_z_close_simunits[stage],2));
        }
    }
    else { //single chromosome
        if (site + 1 == oriC) {
//            if (pol_close_pole(site+1,thread_num)) {
//                energy_change += alpha[stage] * (std::abs(prop_move1[2] - pole[thread_num]) - std::abs(polym[site + 1][2] - pole[thread_num]));
//            }
            energy_change += alpha[stage] * (prop_move1[2] - polym[site + 1][2]);
            energy_change += alpha2[stage] * (pow(prop_move1[2] -xp_z_close_simunits[stage],2) - pow(polym[site + 1][2] -xp_z_close_simunits[stage],2));
        }
    }

    int moved_site1 = (site + 1) % pol_length;

    // energy differences from means: use *distance to nearest pole / cell length*
    if (energy_mean_map[stage].find(moved_site1) != energy_mean_map[stage].end()) {
        double coeff = energy_mean_map[stage][moved_site1];

        double total_cell_length = length[stage] + 2*radius;
        double offset_stage = (int(length[stage]) % 2) / 2.0;
        double left_pole_z  = -(length[stage]/2 + radius - offset_stage);
        double right_pole_z = +(length[stage]/2 + radius - offset_stage);

        auto nearest_norm = [&](double z) {
            double dleft  =  z - left_pole_z;
            double dright =  right_pole_z - z;
            return std::min(dleft, dright) / total_cell_length;   // same metric as normalize_measured_z
        };

        double d_old = nearest_norm(polym[moved_site1][2]);
        double d_new = nearest_norm(prop_move1[2]);

        if (!is_replicated[thread_num][moved_site1]) {
            energy_change += coeff * (d_new - d_old);
        } else if (include_replicated_in_mean) {
            energy_change += coeff * (d_new - d_old) / 2.0;
        }
    }
    //energy differences from separations
    // --- Existing generic separation (ori–ori, etc.) ---
    if (is_replicated[thread_num][moved_site1] &&
        energy_separation_map[stage].find(moved_site1) != energy_separation_map[stage].end())
    {
        energy_change += energy_separation_map[stage][moved_site1] *
            (std::abs(other_polym[moved_site1][2] - prop_move1[2]) -
            std::abs(other_polym[moved_site1][2] - polym[moved_site1][2]));
    }

    // --- NEW: ori–Ter separation at stage 0 (single-chromosome) ---
    if (stage == 0 &&
        energy_separation_map[stage].find(Ter) != energy_separation_map[stage].end())
    {
        double coeff = energy_separation_map[stage][Ter];   // this is energ_coeff_separation[0][1]

        int moved = moved_site1;              // (site+1)%pol_length
        bool moved_is_ori = (moved == oriC);
        bool moved_is_ter = (moved == Ter);

        // Only matters if the move touches oriC or Ter on the ring
        if (!lin_pol && (moved_is_ori || moved_is_ter)) {
            // old positions on the ring
            double z_ori_old = polymer[thread_num][oriC][2];
            double z_ter_old = polymer[thread_num][Ter][2];

            // new positions: start from old, then override the moved one
            double z_ori_new = z_ori_old;
            double z_ter_new = z_ter_old;

            if (moved_is_ori) {
                z_ori_new = prop_move1[2];
            } else if (moved_is_ter) {
                z_ter_new = prop_move1[2];
            }

            double d_old = std::abs(z_ori_old - z_ter_old);
            double d_new = std::abs(z_ori_new - z_ter_new);

            energy_change += coeff * (d_new - d_old);
        }
    }

        return energy_change;
    }

double delta_E_crankshaft(std::vector<Vector3i>& polym, std::vector<Vector3i>& other_polym, int site, Vector3i prop_move1, Vector3i prop_move2, int thread_num, int lin_pol) {
    double energy_change1 = delta_E_other(polym, other_polym, site, prop_move1, thread_num, lin_pol);
    double energy_change2 = delta_E_other(polym, other_polym, site+1, prop_move2, thread_num, lin_pol);
    return energy_change1 + energy_change2;
}

#endif
