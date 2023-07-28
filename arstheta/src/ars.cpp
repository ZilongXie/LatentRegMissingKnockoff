// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp11)]]

#define ARMA_64BIT_WORD 1

#include <Rcpp.h>
#include "logdensity.h"
#include <ctime>
using namespace Rcpp;
using namespace std;


double logsumexp(double u, double v) {
  double rv  = max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)));
  return rv;
}
double logdiffexp(double u, double v) {
  double rv = max(u, v) + log(exp(u - max(u, v)) - exp(v - max(u, v)));
  return rv;
}


struct HullSegment {
public:
  double x;
  double hx;
  double hpx;
  double z_left;
  double z_right;
  double hu_left;
  double hu_right;
  double scum_left;
  double scum_right;
};

class UpperHull {
public:
  int n_segments;
  vector<HullSegment> segments;
  double log_cu;
  double xlb;
  double xrb;
  
  UpperHull() {}; // Default constructor
  
  UpperHull(NumericVector x, NumericVector hx, NumericVector hpx, double _xlb, double _xrb) {
    // x assumed sorted!
    
    // No of segments
    n_segments = x.size();
    
    // Set bounds
    xlb = _xlb;
    xrb = _xrb;
    
    // Ensure left (right) most point is left (right) of mode
    if (hpx.at(0) < 0) {
      stop("Smallest starting point is right of mode.");
    }
    if (hpx.at(hpx.size()-1) > 0) {
      stop("Largest starting point is left of mode.");
    }
    
    // Initialize segments
    for (int i = 0; i < n_segments; ++i) {
      
      HullSegment segment;
      segment.x = x.at(i);
      segment.hx = hx.at(i);
      segment.hpx = hpx.at(i);
      
      segments.push_back(segment);
    }
    
    // Update z and scum
    update_z();
    update_scum();
    
  };
  
  void update_z() {
    
    for (int i = 0; i < n_segments; ++i) {
      
      HullSegment segment = segments.at(i); // We make a copy and reassign to segments below; inefficient, but safer
      
      // Assign left z
      if (i == 0) {
        segment.z_left = xlb;
      } else {
        segment.z_left = segments.at(i-1).z_right;
      }
      
      // Assign right z
      if (i < n_segments - 1) {
        HullSegment next_segment;
        next_segment = segments.at(i+1);
        double num = segment.hx - next_segment.hx + next_segment.hpx*(next_segment.x - segment.x);
        double denom = next_segment.hpx - segment.hpx;
        segment.z_right = segment.x + (num/denom);
      } else {
        segment.z_right = xrb;
      }
      
      // Assign hull values
      segment.hu_left = segment.hx - segment.hpx*(segment.x - segment.z_left);
      segment.hu_right = segment.hx + segment.hpx*(segment.z_right - segment.x);
      
      // Reassign segment copy
      segments.at(i) = segment;
    }
    
  };
  
  void update_scum() {
    
    // Compute normalizer
    vector<double> log_cu_vec(n_segments);
    log_cu = 0;
    for (int i = 0; i < n_segments; ++i) {
      HullSegment segment;
      segment = segments.at(i);
      // double cu_i = (1/segment.hpx)*(exp(segment.hu_right) - exp(segment.hu_left));
      double log_cu_i;
      if (segment.hpx > 0) {
        log_cu_i = -log(segment.hpx) + logdiffexp(segment.hu_right, segment.hu_left);
      } else {
        log_cu_i = logdiffexp(segment.hu_left - log(-segment.hpx), segment.hu_right - log(-segment.hpx));
      }
      
      if (i == 0) {
        log_cu = log_cu_i;
      } else {
        log_cu = logsumexp(log_cu, log_cu_i);
      }
      log_cu_vec.at(i) = log_cu_i;
    }
    
    // Compute and assign scum
    vector<double> scum_vec(n_segments);
    for (int i = 0; i < n_segments; ++i) {
      if (i == 0) {
        scum_vec.at(0) = exp(log_cu_vec.at(0) - log_cu);
        segments.at(0).scum_left = 0;
        segments.at(0).scum_right = scum_vec.at(0);
      } else {
        scum_vec.at(i) = scum_vec.at(i-1) + exp(log_cu_vec.at(i) - log_cu);
        segments.at(i).scum_left = segments.at(i-1).scum_right;
        segments.at(i).scum_right = scum_vec.at(i);
      }
    }
  };
  
  double sample() {
    
    double x_sample;
    double u = R::unif_rand();
    double min_scum = segments.at(0).scum_right;
    if (u < min_scum) {
      // Sample is on the left-most segment
      HullSegment s = segments.at(0);
      
      // Draw a sample from the upper hull
      
      // x_sample = s.z_right + (1/s.hpx)*log(1 + (s.hpx*exp(log_cu)*(u-s.scum_right))/exp(s.hu_right));
      // ^ This often leads to over/underflow
      
      if (s.hpx*(u-s.scum_right) > 0) { // LH term might be negative, so we can't compute the log naively
        x_sample = s.z_right + (1/s.hpx)*logsumexp(0, log(s.hpx*(u-s.scum_right)) + log_cu - s.hu_right);
      } else {
        x_sample = s.z_right + (1/s.hpx)*logdiffexp(0, log(-1*s.hpx*(u - s.scum_right)) + log_cu - s.hu_right);
      }
      
    } else {
      // Determine which segment the sample is on
      int segment_id;
      for (int i = 1; i < n_segments; ++i) {
        if (u > segments.at(i).scum_left && u < segments.at(i).scum_right) {
          segment_id = i;
        }
      }
      HullSegment s = segments.at(segment_id);
      
      // Draw a sample from the upper hull
      
      // x_sample = s.z_left + (1/s.hpx)*log(1 + (s.hpx*exp(log_cu)*(u-s.scum_left))/exp(s.hu_left));
      // ^ This is prone to over/underflow
      
      if (s.hpx > 0) { // s.hpx might be negative, so we can't compute the log naively!
        x_sample = s.z_left + (1/s.hpx)*logsumexp(0, log(s.hpx) + log_cu + log(u - s.scum_left) - s.hu_left);
      } else {
        x_sample = s.z_left + (1/s.hpx)*logdiffexp(0, log(-s.hpx) + log_cu + log(u - s.scum_left) - s.hu_left);
      }
    }
    
    return x_sample;
  };
  
  void add_segment(double x, double hx, double hpx) {
    
    HullSegment segment;
    segment.x = x;
    segment.hx = hx;
    segment.hpx = hpx;
    
    // Determine segment position
    int iter = 0;
    while(iter < n_segments) {
      if (x < segments.at(iter).x) {
        break;
      }
      iter++;
    }
    
    // Insert segment
    segments.insert(segments.begin()+iter, segment);
    n_segments = segments.size();
    
    // Update z and scum
    update_z();
    update_scum();
    
  };
  
  double get_hu(double x) {
    
    // Determine which segment x lies on
    int segment_id;
    for (int i = 0; i < n_segments; ++i) {
      if (x > segments.at(i).z_left && x <= segments.at(i).z_right) {
        segment_id = i;
        break;
      }
    }
    HullSegment s = segments.at(segment_id);
    
    // Get hu
    double hu;
    if (segment_id == 0) {
      hu = s.hu_right - (s.z_right - x)*s.hpx;
    } else {
      hu = s.hu_left + (x - s.z_left)*s.hpx;
    }
    
    return hu;
  };
};

class LowerHull {
public:
  vector<double> x;
  vector<double> hx;
  vector<double> hpx;
  int n_points;
  
  LowerHull() {}; // Default constructor
  
  LowerHull(NumericVector _x, NumericVector _hx, NumericVector _hpx) {
    // x assumed sorted!
    
    x = as<vector<double> >(_x);
    hx = as<vector<double> >(_hx);
    hpx = as<vector<double> >(_hpx);
    n_points = x.size();
    
  };
  
  void add_segment(double new_x, double new_hx, double new_hpx) {
    
    // Determine point position
    int iter = 0;
    while (iter < n_points) {
      if (new_x < x.at(iter)) {
        break;
      }
      iter++;
    }
    
    // Assign
    x.insert(x.begin()+iter, new_x);
    hx.insert(hx.begin()+iter, new_hx);
    hpx.insert(hpx.begin()+iter, new_hpx);
    
    n_points = x.size();
  }
  
  double get_hl(double _x) {
    
    // Determine point position
    int iter = 0;
    while (iter < n_points) {
      if (_x < x.at(iter)) {
        break;
      }
      iter++;
    }
    
    double rv;
    if (iter == 0) {
      rv = -numeric_limits<double>::infinity();
    } else if (iter == n_points) {
      rv = -numeric_limits<double>::infinity();
    } else {
      double x_left = x.at(iter-1);
      double x_right = x.at(iter);
      double hx_left = hx.at(iter-1);
      double hx_right = hx.at(iter);
      double d = (hx_right-hx_left)/(x_right-x_left);
      rv = hx_left + (_x-x_left)*d;
    }
    
    return rv;
  }
  
};

class ARS {
  // AR sampling using an object of class LogDensity
public:
  UpperHull upper_hull;
  LowerHull lower_hull;
  double xlb;
  double xrb;
  int max_points;
  
  const static int MAX_REJECTIONS = 100000;
  
  ARS() {}; // Default constructor
  
  ARS(LogDensity* const log_density, NumericVector x, double xlb, double xrb, int max_points) :
    xlb(xlb),
    xrb(xrb),
    max_points(max_points) {
    
    // Check bounds
    if (xlb >= xrb) {
      stop("Upper bound is not largetr than lower bound.");
    }
    
    // We need at least two starting points
    if (x.size() < 2) {
      stop("At least two starting points required.");
    }
    
    // Order x
    x.sort();
    
    // Check points in bounds
    if (x.at(0) < xlb || x.at(x.size()-1) > xrb) {
      stop("Starting point out of bound.");
    }
    
    // Evaluate funcs
    NumericVector hx = log_density->h(x);
    NumericVector hpx = log_density->h_prime(x);
    
    // Try to make starting points valid (if they're invalid)
    int max_tries = 10;
    if (hpx.at(0) <= 0) {
      // Move left-most point left until valid
      
      double hpx_left = hpx.at(0);
      double x_left = x.at(0);
      int tries = 0;
      while (hpx_left <= 0 && tries < max_tries) {
        
        if (isfinite(xlb)) {
          x_left -= (x_left-xlb)/2; // Move half-way to the limit
        } else {
          x_left -= pow(2, tries); // Move left by power of 2
        }
        
        hpx_left = log_density->h_prime(x_left);
        
        tries++;
      }
      
      if (tries < max_tries) {
        hpx.at(0) = hpx_left;
        x.at(0) = x_left;
        hx.at(0) = log_density->h(x_left);
      } else {
        stop("Could not find valid lower starting point.");
      }
    }
    
    int last_ind = hpx.size() - 1;
    if (hpx.at(last_ind) >= 0) {
      // Move right-most point right until valid
      
      double hpx_right = hpx.at(last_ind);
      double x_right = x.at(last_ind);
      int tries = 0;
      while (hpx_right >= 0 && tries < max_tries) {
        
        if (isfinite(xrb)) {
          x_right += (xrb - x_right)/2; // Move half-way to the limit
        } else {
          x_right += pow(2, tries); // Move right by power of 2
        }
        
        hpx_right = log_density->h_prime(x_right);
        
        tries++;
      }
      
      if (tries < max_tries) {
        hpx.at(last_ind) = hpx_right;
        x.at(last_ind) = x_right;
        hx.at(last_ind) = log_density->h(x_right);
      } else {
        Rcout << "x candidates: " << x << endl;
        stop("Could not find valid upper starting point.");
      }
    }
    
    // Create the hull
    upper_hull = UpperHull(x, hx, hpx, xlb, xrb);
    lower_hull = LowerHull(x, hx, hpx);
  };
  
  NumericVector sample(unsigned N, LogDensity* const log_density) {
    
    vector<double> samples;
    int rejections = 0;
    while(samples.size() < N) {
      
      double x_sample = upper_hull.sample();
      double u = R::runif(0, 1);
      if (u < exp(lower_hull.get_hl(x_sample) - upper_hull.get_hu(x_sample))) {
        // Accept!
        samples.push_back(x_sample);
        rejections = 0;
      } else {
        
        double hx = log_density->h(x_sample);
        
        if (u < exp(hx - upper_hull.get_hu(x_sample))) {
          // Accept!
          samples.push_back(x_sample);
          rejections = 0;
        } else {
          // Reject!
          rejections++;
        }
        
        // Add hull segment
        int points = lower_hull.x.size();
        if (points < max_points) {
          double hpx = log_density->h_prime(x_sample);
          upper_hull.add_segment(x_sample, hx, hpx);
          lower_hull.add_segment(x_sample, hx, hpx);
        }
      }
      
      if (rejections > MAX_REJECTIONS) {
        // Rcout << "Warning: Maximum number of rejections reached. Returning current sample." << endl;
        // samples.push_back(x_sample);
        Rcout << "Error: Maximum number of rejections reached." << endl;
        stop("Maximum number of rejections reached");
      }
      
      checkUserInterrupt();
    }
    
    return wrap(samples);
  };
};

class ARS_Rwrapper {
  // A wrapper class around ARS for use from R
public:
  ARS ars;
  RLogDensity r_log_density;
  
  ARS_Rwrapper(Function h, Function h_prime, NumericVector x, double xlb, double xrb, int max_points) :
    r_log_density(h, h_prime) {
    ARS new_ars(&r_log_density, x, xlb, xrb, max_points);
    ars = new_ars;
  };
  
  NumericVector sample(int N) {
    return ars.sample(N, &r_log_density);
  };
  
  NumericVector get_x() {
    return wrap(ars.lower_hull.x);
  };
  
  double get_hu(double x) {
    return ars.upper_hull.get_hu(x);
  };
  
  double get_hl(double x) {
    return ars.lower_hull.get_hl(x);
  };
  
  double get_cu() {
    return exp(ars.upper_hull.log_cu);
  }
};


class ThetaPosterior{
  // Class for sampling from posterior of theta given Y and X
  // p(theta_i | Y_i, X_i) \propto p(Y_i | theta_i) * p(theta_i | X_i)
public:
  NumericMatrix Y;
  NumericVector avec;
  NumericVector d1vec;
  NumericVector d2vec;
  NumericVector mu;
  double var;
  
  vector<ThetaLogPost> densities;
  int N;
  int J;
  
  ThetaPosterior(NumericMatrix Y, NumericVector avec, NumericVector d1vec, NumericVector d2vec, NumericVector mu, double var):
    Y(Y),
    avec(avec),
    d1vec(d1vec),
    d2vec(d2vec),
    mu(mu),
    var(var){
    
    N = Y.nrow();
    J = Y.ncol();
    
    // Check arguments
    if (avec.size() != J || d1vec.size() != J || d2vec.size() != J) {
      stop("Y_i, avec, d1vec and d2vec must have same length.");
    }
    if (mu.size() != N) {
      stop("length of mu does not equal to sample size.");
    }
    if (var <= 1e-4) {
      stop("Input variance is too small!");
    }
    
    // Make the densities
    for (int i = 0; i < N; ++i) {
      // log posterior density of theta
      ThetaLogPost log_density(Y(i, _), avec, d1vec, d2vec, mu.at(i), var);
      densities.push_back(log_density);
    }
  };
  
  NumericVector sample(NumericVector theta_init) {
    
    // Make the samplers
    vector<ARS> samplers;
    for (int i = 0; i < N; ++i) {
      // Initialize sampler with 5 starting points, with at most 100 points in hulls.
      NumericVector theta_init_i = NumericVector::create(theta_init.at(i) - 2, 
                                                         theta_init.at(i),
                                                         theta_init.at(i) + 1,
                                                         theta_init.at(i) + 3,
                                                         theta_init.at(i) + 5);
      ARS ars(&densities.at(i), theta_init_i, theta_init.at(i)-100, theta_init.at(i)+100, 100);
      samplers.push_back(ars);
    };
    
    NumericVector out(N);
    for (int i = 0; i < N; ++i) {
      // only one sample
      NumericVector smp_i = samplers.at(i).sample(1, &densities.at(i));
      out.at(i) = smp_i.at(0);
    }
    return out;
  };
};


RCPP_MODULE(hull) {
  class_<ARS_Rwrapper>("ARS_Rwrapper")
  .constructor<Function, Function, NumericVector, double, double, int>()
  .method("sample", &ARS_Rwrapper::sample)
  .method("get_x", &ARS_Rwrapper::get_x)
  .method("get_hu", &ARS_Rwrapper::get_hu)
  .method("get_hl", &ARS_Rwrapper::get_hl)
  .method("get_cu", &ARS_Rwrapper::get_cu)
  ;
  
  class_<ThetaPosterior>("ThetaPosterior")
    .constructor<NumericMatrix, NumericVector, NumericVector, NumericVector, NumericVector, double>()
    .method("sample", &ThetaPosterior::sample)
  ;
}
