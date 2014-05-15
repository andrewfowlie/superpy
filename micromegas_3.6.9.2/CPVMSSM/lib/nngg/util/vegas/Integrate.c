/*
	Integrate.c
		integrate over the unit hypercube
		this file is part of Vegas
		last modified 16 Jul 04
*/


static int Integrate(creal epsrel, creal epsabs,
  cint flags, ccount mineval, ccount maxeval,
  ccount nstart, ccount nincrease,
  real *integral, real *error, real *prob)
{
  Cumulants cumul[NCOMP];
  Grid *grid;
  real *sample;
  count dim, comp, bin, df;
  int fail = 1;
  count nsamples = nstart;

  if( VERBOSE > 1 ) {
    char s[256];
    sprintf(s, "Vegas input parameters:\n"
      "  epsrel %g\n  epsabs %g\n"
      "  flags %d\n  mineval %d\n  maxeval %d\n"
      "  nstart %d\n  nincrease %d",
      epsrel, epsabs, flags, mineval, maxeval, nstart, nincrease);
    Print(s);
  }

#ifdef MLVERSION
  if( setjmp(abort_) ) goto abort;
#endif

  IniRandom(2*maxeval, ndim_);
  grid = GetGrid();
  Zap(cumul);

  /* main iteration loop */

  for( df = 0; ; ++df ) {
    count n;
    creal jacobian = 1./nsamples;
    real *w, *x, *f, *lastf;
    bin_t *bins;

    SamplesAlloc(sample, nsamples);
    w = sample;
    x = w + nsamples;
    f = x + nsamples*ndim_;
    lastf = f + nsamples*ncomp_;
    bins = (bin_t *)lastf;

    for( n = nsamples; n; --n ) {
      real weight = jacobian;

      GetRandom(x);

      for( dim = 0; dim < ndim_; ++dim ) {
        creal pos = *x*NBINS;
        ccount bin = pos;
        creal prev = (bin == 0) ? 0 : grid[dim][bin - 1];
        creal diff = grid[dim][bin] - prev; 
        *x++ = prev + (pos - bin)*diff;
        *bins++ = bin;
        weight *= diff*NBINS;
      }

      *w++ = weight;
    }

    DoSample(nsamples, w, f);

    w = sample;

    while( f < lastf ) {
      creal weight = *w++;

      for( comp = 0; comp < ncomp_; ++comp ) {
        Cumulants *c = &cumul[comp];
        creal wfun = weight*(*f++);
        c->sum += wfun;
        c->sqsum += Sq(wfun);
      }
    }

    fail = 0;

    /* compute the integral and error values */

    for( comp = 0; comp < ncomp_; ++comp ) {
      Cumulants *c = &cumul[comp];
      real avg, sigsq;
      real w = Weight(c->sum, c->sqsum, nsamples);

      sigsq = 1/(c->weightsum += w);
      avg = sigsq*(c->avgsum += w*c->sum);

      c->avg = LAST ? (sigsq = 1/w, c->sum) : avg;
      c->err = sqrt(sigsq);
      fail |= (c->err > Max(epsabs, fabs(c->avg)*epsrel));

      if( df == 0 ) c->guess = c->sum;
      else {
        c->chisum += w *= c->sum - c->guess;
        c->chisqsum += w*c->sum;
      }
      c->chisq = c->chisqsum - avg*c->chisum;

      c->sum = c->sqsum = 0;
    }

    if( VERBOSE ) {
      char s[128 + 128*NCOMP], *p = s;

      p += sprintf(p, "\n"
        "Iteration %d:  %d integrand evaluations so far",
        df + 1, neval_);

      for( comp = 0; comp < ncomp_; ++comp ) {
        cCumulants *c = &cumul[comp];
        p += sprintf(p, "\n[%d] %g +- %g  \tchisq %g (%d df)",
          comp + 1, c->avg, c->err, c->chisq, df);
      }

      Print(s);
    }

    if( fail == 0 && neval_ >= mineval ) break;

    if( neval_ >= maxeval ) break;

    Reweight(grid, sample, x, f, cumul);

    free(sample);
    nsamples += nincrease;
  }

  for( comp = 0; comp < ncomp_; ++comp ) {
    cCumulants *c = &cumul[comp];
    integral[comp] = c->avg;
    error[comp] = c->err;
    prob[comp] = ChiSquare(c->chisq, df);
  }

abort:
  free(sample);
  PutGrid(grid);

  return fail;
}
