#pragma once

#include <random>
#include <utility>
#include <vector>
#include <stack>
#include <iostream>

// Implementation of Vose's alias method, given a list of probabilities
// representing a finite distribution builds a datastructure in O(n) time
// that allows O(1) drawing from the distribution.
template<size_t N>
class DiscreteNDistribution
{
private:
    std::array<double, N> probabilities_;
    std::array<int, N> alias_;

    // std::uniform_real_distribution<double> biased_coin_;
    // std::uniform_int_distribution<int> fair_die_;


public:
    DiscreteNDistribution<N>(const std::vector<double> &probabilities, double normalizingFactor=1.0) {
        assert(probabilities.size() == N);  // safety check

        std::stack<std::pair<int, double>> small_;
        std::stack<std::pair<int, double>> large_;


        for(int i = 0; i < N; i++) {
            double scaled_prob = N*probabilities[i]*normalizingFactor;
            std::pair<int, double> current_prob(i, scaled_prob);
            if (current_prob.second < 1) {
                small_.push(current_prob);
            } else {
                large_.push(current_prob);
            }
        }
        while (!small_.empty() && !large_.empty()) {
            std::pair<int, double> s = small_.top();
            std::pair<int, double> l = large_.top();
            small_.pop();
            large_.pop();

            probabilities_[s.first] = s.second;
            alias_[s.first] = l.first;

            l.second = (l.second + s.second) - 1;
            if(l.second < 1) {
                small_.push(l);
            } else {
                large_.push(l);
            }
        }
        while (!large_.empty()) {
            std::pair<int, double> l = large_.top();
            large_.pop();
            probabilities_[l.first] = 1;
        }
        while (!small_.empty()) {
            std::pair<int, double> s = small_.top();
            small_.pop();
            probabilities_[s.first] = 1;
        }
    }

    template<typename RngType = std::mt19937_64>
    int drawSample(RngType &rng) {
        uint64_t random_bits = rng();
    int die_roll = random_bits % N;  // Uses lower ~5 bits for N=20
    double coin_flip = ((random_bits >> 8) & 0xFFFFFFFF) * 2.3283064365386963e-10;  // Uses different bits
        
        if (coin_flip < probabilities_[die_roll]) return die_roll + 1;
        return alias_[die_roll] + 1;
    }

    void printTable() {
        for(auto &i: probabilities_){
            std::cout << i << " ";
        }
        std::cout << "\n";
        for(auto &i: alias_){
            std::cout << i << " ";
        }
        std::cout << "\n";
    }

    // static void setSeed(int seed) {
    //     rng_ = std::mt19937_64(seed);
    // }

    std::vector<std::pair<double, int>> getTable() {
        std::vector<std::pair<double, int>> prob_alias_table;
        for(int i=0; i < probabilities_.size(); i++){
            std::pair<double, int> current_item (probabilities_[i], alias_[i]);
            prob_alias_table.push_back(current_item);
        }
        return prob_alias_table;
    }

    double getProb(int i) {
        return probabilities_[i];
    }
};