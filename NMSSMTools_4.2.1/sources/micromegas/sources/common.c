#include "micromegas.h"

MOcommonSTR mocommon_=
{{
  0        , /* 0      null         */
  100      , /* 1      Mcdm         */
  0.0191   , /* 2      ScalarFFPd   */
  0.0153   , /* 3      ScalarFFPu   */ 
  0.0447   , /* 4      ScalarFFPs   */
  -0.427   , /* 5      pVectorFFPd  */
   0.842   , /* 6      pVectorFFPu  */
  -0.085   , /* 7      pVectorFFPs  */
  -0.23    , /* 8      SigmaFFPd    */
   0.84    , /* 9      SigmaFFPu    */
   -0.046  , /* 10     SigmaFFPs    */

  0.0273   , /* 11     ScalarFFNd   */
  0.0110   , /* 12     ScalarFFNu   */
  0.0447   , /* 13     ScalarFFNs   */
  0.842    , /* 14     pVectorFFNd  */
  -0.427   , /* 15     pVectorFFNu  */
  -0.085   , /* 16     pVectorFFNs  */
  0.84     , /* 17     SigmaFFNd    */
  -0.23    , /* 18     SigmaFFNu    */
  -0.046   , /* 19     SigmaFFNs    */
   0.52    , /* 20     Fermi_a      */
   -0.6    , /* 21     Fermi_b      */
   1.23    , /* 22     Fermi_c      */
   8.0     , /* 23     Rsun         */
   0.3     , /* 24     rhoDM(Sun)   */
  225.2    , /* 25     Vearth       */ 
   0.0112  , /* 26     K_dif        */
   4       , /* 27     L_dif        */
   0.7     , /* 28     Delta_dif    */
   1.E16   , /* 29     Tau_dif      */
   12      , /* 30     Vc_dif       */
   20      , /* 31     Rdisk        */
   0       , /* 32     deltaY  (abandence asymmtry) */
   0       , /* 33     dmAsymm (log(dm/dm_bar)) */
   600     , /* 34     Vesc         */
   220      /* 35     Vrot         */
//   225       /* 36     Vearth       */
}};

void forceug_(int * key) { ForceUG=*key;}
