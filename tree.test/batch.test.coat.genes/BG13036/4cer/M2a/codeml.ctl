      seqfile = cds2.phylip   * sequence data filename
     treefile = tree.cereus.M2a.nwk      * tree structure file name
      outfile = results.cereus.M2a.txt   * main result file name

        noisy = 9      * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1      * 1:detailed output
      runmode = 0      * 0:user defined tree

      seqtype = 1      * 1:codons
    CodonFreq = 2      * 0:equal, 1:F1X4, 2:F3X4, 3:F61

       model = 0    
                       * 1:separate omega for each branch
                       * 2:user specified dN/dS ratios for branches

      NSsites = 2     * 
                    * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;
                    * 4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal

        icode = 0      * 0:universal code

    fix_kappa = 0      * 1:kappa fixed, 0:kappa to be estimated
        kappa = 2      * initial or fixed kappa

    fix_omega = 0      * 1:omega fixed, 0:omega to be estimated 
        omega = 0.2    * initial omega

                       * comments:
                       * H0 in Table 3: model = 0
                       * H1 in Table 3: model = 2
                       * H2 in Table 3: model = 2
                       * H3 in Table 3: model = 2
