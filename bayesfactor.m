function [bayesfactor,likelihoodtheory,likelihoodnull] = bayesfactor(samplemean, sampleSE, uniform, varargin)
%bayesfactor: compute bayes factor on difference of means
%
%   [bayesfactor] = bayesfactor(samplemean,sampleSE,1,lower,upper) computes
%   the bayes factor that the samplemean and sampleSE differences between
%   conditions lie within a uniform distribution between the lower and upper bounds of plausible values
%   for the alternative hypothesis
%   The lower limit can often be set to zero, but when there is a minimal value below which the effect is too small for the theory to be either true or interesting, that can form the lower limit. If there is no theoretical constraint on the upper limit, but the scale has a natural upper limit (e.g. it is a rating scale from 0-7, then a difference between conditions cannot exceed 7), then that can be specified as the upper limit, if you really have no grounds for otherwise limiting it.
%
%   [bayesfactor] =
%   bayesfactor(samplemean,sampleSE,0,theorymean,theorySE,tail) computes
%   the bayes factor that the samplemean and sampleSE differences between
%   conditions lie within a normal distribution with theorymean and
%   theorySE as predicted by the alternative hypothesis.
%
%   [bayesfactor,likelihoodtheory,likelihoodnull] = bayesfactor(...) also
%   returns the likelihoods for both theories.
%
%   A Bayes factor of 3 or more can be taken as substantial evidence for
%   your theory (and against the null).
%   A Bayes factor of 1/3 or less can be taken as evidence for the null (and against your theory).
%   Bayes factors between 1/3 and 3 show the data do not provide much evidence to distinguish your theory from the null.
%
%   Algorithm and documentation from:
%   http://www.lifesci.sussex.ac.uk/home/Zoltan_Dienes/inference/Bayes.htm


    sd = sampleSE; 
    sd2 = sd*sd; 
    obtained = samplemean; 

    %uniform = input('is the distribution of p(population value|theory) uniform? 1= yes 0=no '); 

    if uniform
        lower = varargin{1};%input('What is the lower bound? '); 
        upper = varargin{2};%('What is the upper bound? '); 
    else
        meanoftheory = varargin{1};%input('What is the mean of p(population value|theory)? '); 
        sdtheory = varargin{2};%input('What is the standard deviation of p(population value|theory)? '); 
        omega = sdtheory*sdtheory;    
        tail = varargin{3};% input('is the distribution one-tailed or two-tailed? (1/2) '); 
    end

    area = 0; 
    if uniform
        theta = lower; 
        incr = (upper- lower)/2000; 
    else
        theta = meanoftheory - 5*(omega)^0.5;
        incr =  (omega)^0.5/200;
    end

    %define the normal function
    %normaly = @(mn, variance, x) 2.718283^(- (x - mn)*(x - mn)/(2*variance))/realsqrt(2*pi*variance); 
    normaly = @(mn, variance, x) exp(1)^(- (x - mn)*(x - mn)/(2*variance))/realsqrt(2*pi*variance); 

    for i = -1000:1000
        theta = theta + incr; 
        if uniform
            dist_theta = 0; 
            if and(theta >= lower, theta <= upper) 
                dist_theta = 1/(upper-lower); 
            end               
        else %distribution is normal 
            if tail == 2 
                dist_theta = normaly(meanoftheory, omega, theta); 
            else 
                dist_theta = 0; 
                if theta > 0 
                    dist_theta = 2*normaly(meanoftheory, omega, theta); 
                end 
            end 
        end 
        height = dist_theta * normaly(theta, sd2, obtained); %p(population value=theta|theory)*p(data|theta) 
        area = area + height*incr; %integrating the above over theta 
    end 

    likelihoodtheory = area;
    likelihoodnull = normaly(0, sd2, obtained);
    bayesfactor = likelihoodtheory/likelihoodnull;
end