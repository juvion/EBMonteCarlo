double MonteCarloClass::SimpleSpringEnergy(RNA *rna1, RNA *rna2, RNA *rna3) {
    int dist_a = 0;
    int dist_b = 0;
    //iterate the each bases of the sequence.
    for (int i = 1; i <= rna1->GetSequenceLength(); i++) {
        //find the pair to base index for rna1, rna2 and rna3
        int pair_to1 = rna1->GetPair(i);
        int pair_to2 = rna2->GetPair(i);
        int pair_to3 = rna3->GetPair(i);
        //measure the base pair distance between rna1 and rna2
        if (!(pair_to1 == 0 && pair_to2 ==0)) {
            if (pair_to1 != pair_to2) {
                if (pair_to1 == 0 || pair_to2 == 0) {
                    dist_a++;
                }
                else {
                    dist_a ++;
                    dist_a ++; 
                }
            }
        }
        //measure the base pair distance between rna2 and rna3
        if (!(pair_to2 == 0 && pair_to3 ==0)) {
            if (pair_to2 != pair_to3) {
                if (pair_to2 == 0 || pair_to3 == 0) {
                    dist_b++;
                }
                else {
                    dist_b ++;
                    dist_b ++;
                }
            }
        }
    }
    //iterate the each bases of the sequence.
    //define the 
    double Spring_energy = SpringConst * (double(dist_a * dist_a) + double(dist_b * dist_b));
    return Spring_energy;
}