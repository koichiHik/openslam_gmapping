#ifndef GRIDSLAMPROCESSOR_H
#define GRIDSLAMPROCESSOR_H

#include <climits>
#include <limits>
#include <fstream>
#include <vector>
#include <deque>
#include <gmapping/particlefilter/particlefilter.h>
#include <gmapping/utils/point.h>
#include <gmapping/utils/macro_params.h>
#include <gmapping/log/sensorlog.h>
#include <gmapping/sensor/sensor_range/rangesensor.h>
#include <gmapping/sensor/sensor_range/rangereading.h>
#include <gmapping/scanmatcher/scanmatcher.h>
#include "motionmodel.h"


namespace GMapping {

  /**This class defines the basic GridFastSLAM algorithm.  It
     implements a rao blackwellized particle filter. Each particle
     has its own map and robot pose.<br> This implementation works
     as follows: each time a new pair odometry/laser reading is
     received, the particle's robot pose is updated according to the
     motion model.  This pose is subsequently used for initalizing a
     scan matching algorithm.  The scanmatcher performs a local
     optimization for each particle.  It is initialized with the
     pose drawn from the motion model, and the pose is corrected
     according to the each particle map.<br>
     In order to avoid unnecessary computation the filter state is updated 
     only when the robot moves more than a given threshold.
  */
  class GridSlamProcessor{
  public:

    
    /**This class defines the the node of reversed tree in which the trajectories are stored.
       Each node of a tree has a pointer to its parent and a counter indicating the number of childs of a node.
       The tree is updated in a way consistent with the operation performed on the particles.
   */
    struct TNode{
      /**Constructs a node of the trajectory tree.
       @param pose:      the pose of the robot in the trajectory
       @param weight:    the weight of the particle at that point in the trajectory
       @param accWeight: the cumulative weight of the particle
       @param parent:    the parent node in the tree
       @param childs:    the number of childs
      */
      TNode(const OrientedPoint& pose, double weight, TNode* parent=0, unsigned int childs=0);

      /**Destroys a tree node, and consistently updates the tree. If a node whose parent has only one child is deleted,
       also the parent node is deleted. This because the parent will not be reacheable anymore in the trajectory tree.*/
      ~TNode();

      /**The pose of the robot*/
      OrientedPoint pose; 
      
      /**The weight of the particle*/
      double weight;

      /**The sum of all the particle weights in the previous part of the trajectory*/
      double accWeight;

      double gweight;


      /**The parent*/
      TNode* parent;

      /**The range reading to which this node is associated*/
      const RangeReading* reading;

      /**The number of childs*/
      unsigned int childs;

      /**counter in visiting the node (internally used)*/
      mutable unsigned int visitCounter;

      /**visit flag (internally used)*/
      mutable bool flag;
    };
    
    typedef std::vector<GridSlamProcessor::TNode*> TNodeVector;
    typedef std::deque<GridSlamProcessor::TNode*> TNodeDeque;
    
    /**This class defines a particle of the filter. Each particle has a map, a pose, a weight and retains the current node in the trajectory tree*/
    struct Particle{
      /**constructs a particle, given a map
	 @param map: the particle map
      */
      Particle(const ScanMatcherMap& map);

      /** @returns the weight of a particle */
      inline operator double() const {return weight;}
      /** @returns the pose of a particle */
      inline operator OrientedPoint() const {return pose;}
      /** sets the weight of a particle
	  @param w the weight
      */
      inline void setWeight(double w) {weight=w;}
      /** The map */
      ScanMatcherMap map;
      /** The pose of the robot */
      OrientedPoint pose;

      /** The pose of the robot at the previous time frame (used for computing thr odometry displacements) */
      OrientedPoint previousPose;

      /** The weight of the particle */
      double weight;

      /** The cumulative weight of the particle */
      double weightSum;

      double gweight;

      /** The index of the previous particle in the trajectory tree */
      int previousIndex;

      /** Entry to the trajectory tree */
      TNode* node; 
    };
	
    
    typedef std::vector<Particle> ParticleVector;
    
    /** Constructs a GridSlamProcessor, initialized with the default parameters */
    GridSlamProcessor();

    /** Constructs a GridSlamProcessor, whose output is routed to a stream.
     @param infoStr: the output stream
    */
    GridSlamProcessor(std::ostream& infoStr);
    
    /** @returns  a deep copy of the grid slam processor with all the internal structures.
    */
    GridSlamProcessor* clone() const;
    
    /**Deleted the gridslamprocessor*/
    virtual ~GridSlamProcessor();
    
    //methods for accessing the parameters
    void setSensorMap(const SensorMap& smap);
    void init(unsigned int size, double xmin, double ymin, double xmax, double ymax, double delta, 
	      OrientedPoint initialPose=OrientedPoint(0,0,0));
    void setMatchingParameters(double urange, double range, double sigma, int kernsize, double lopt, double aopt, 
			       int iterations, double likelihoodSigma=1, double likelihoodGain=1, unsigned int likelihoodSkip=0);
    void setMotionModelParameters(double srr, double srt, double str, double stt);
    void setUpdateDistances(double linear, double angular, double resampleThreshold);
    void setUpdatePeriod(double p) {period_=p;}
    
    //the "core" algorithm
    void processTruePos(const OdometryReading& odometry);
    bool processScan(const RangeReading & reading, int adaptParticles=0);
    
    /**This method copies the state of the filter in a tree.
     The tree is represented through reversed pointers (each node has a pointer to its parent).
     The leafs are stored in a vector, whose size is the same as the number of particles.
     @returns the leafs of the tree
    */
    TNodeVector getTrajectories() const;
    void integrateScanSequence(TNode* node);
    
    /**the scanmatcher algorithm*/
    ScanMatcher m_matcher;
    /**the stream used for writing the output of the algorithm*/
    std::ofstream& outputStream();
    /**the stream used for writing the info/debug messages*/
    std::ostream& infoStream();
    /**@returns the particles*/
    inline const ParticleVector& getParticles() const {return m_particles; }
    
    inline const std::vector<unsigned int>& getIndexes() const{return m_indexes; }
    int getBestParticleIndex() const;
    //callbacks
    virtual void onOdometryUpdate();
    virtual void onResampleUpdate();
    virtual void onScanmatchUpdate();

  //accessor methods
  public:
    /**the maxrange of the laser to consider */
    double getlaserMaxRange() const { return m_matcher.getlaserMaxRange(); }
    void setlaserMaxRange(double laserMaxRange) { m_matcher.setlaserMaxRange(laserMaxRange); }

    /**the maximum usable range of the laser. A beam is cropped to this value. [scanmatcher]*/
    double getusableRange() const { return m_matcher.getusableRange(); }
    void setusableRange(double usableRange) { m_matcher.setusableRange(usableRange); }

    /**The sigma used by the greedy endpoint matching. [scanmatcher]*/
    double getgaussianSigma() const { m_matcher.getgaussianSigma(); }
    void setgaussianSigma(double gaussianSigma) { m_matcher.setgaussianSigma(gaussianSigma); }

    /**The sigma  of a beam used for likelihood computation [scanmatcher]*/
    double getlikelihoodSigma() const { return m_matcher.getlikelihoodSigma(); }
    void setlikelihoodSigma(double likelihoodSigma) {m_matcher.setlikelihoodSigma(likelihoodSigma); }

    /**The kernel in which to look for a correspondence[scanmatcher]*/
    int getkernelSize() const { return m_matcher.getkernelSize(); }
    void setkernelSize(int kernelSize) { m_matcher.setkernelSize(kernelSize); }

    /**The optimization step in rotation [scanmatcher]*/
    double getoptAngularDelta() const { return m_matcher.getoptAngularDelta(); }
    void setoptAngularDelta(double optAngularDelta) { m_matcher.setoptAngularDelta(optAngularDelta); }

    /**The optimization step in translation [scanmatcher]*/
    double getoptLinearDelta() const { return m_matcher.getoptLinearDelta(); }
    void setoptLinearDelta(double optLinearDelta) {m_matcher.setoptLinearDelta(optLinearDelta); }

    /**The number of iterations of the scanmatcher [scanmatcher]*/
    unsigned int getoptRecursiveIterations() const { return m_matcher.getoptRecursiveIterations(); }
    void setoptRecursiveIterations(unsigned int optRecursiveIterations) { m_matcher.setoptRecursiveIterations(optRecursiveIterations); }

    /**the beams to skip for computing the likelihood (consider a beam every likelihoodSkip) [scanmatcher]*/
    unsigned int getlikelihoodSkip() const { return m_matcher.getlikelihoodSkip(); }
    void setlikelihoodSkip(unsigned int likelihoodSkip) { m_matcher.setlikelihoodSkip(likelihoodSkip); }

    /**translational sampling range for the likelihood [scanmatcher]*/
    double getllsamplerange() const { return m_matcher.getllsamplerange(); }
    void setllsamplerange(double llsamplerange) { m_matcher.setllsamplerange(llsamplerange); }

    /**angular sampling range for the likelihood [scanmatcher]*/
    double getlasamplerange() const { return m_matcher.getlasamplerange(); }
    void setlasamplerange(double lasamplerange) { m_matcher.setlasamplerange(lasamplerange); }

    /**translational sampling range for the likelihood [scanmatcher]*/
    double getllsamplestep() const { return m_matcher.getllsamplestep(); }
    void setllsamplestep(double llsamplestep) { m_matcher.setllsamplestep(llsamplestep); }

    /**angular sampling step for the likelihood [scanmatcher]*/
    double getlasamplestep() const {return m_matcher.getlasamplestep(); }
    void setlasamplestep(double lasamplestep) { m_matcher.setlasamplestep(lasamplestep); }

    /**generate an occupancy grid map [scanmatcher]*/
    bool getgenerateMap() const { return m_matcher.getgenerateMap(); }
    void setgenerateMap(bool generateMap) { m_matcher.setgenerateMap(generateMap); }

    /**enlarge the map when the robot goes out of the boundaries [scanmatcher]*/
    bool getenlargeStep() const { return m_matcher.getenlargeStep(); }
    void setenlargeStep(bool enlargeStep) { m_matcher.setenlargeStep(enlargeStep); }

    /**pose of the laser wrt the robot [scanmatcher]*/
    OrientedPoint getlaserPose() const { return m_matcher.getlaserPose(); }
    void setlaserPose(OrientedPoint laserPose) { m_matcher.setlaserPose(laserPose); }

    /**odometry error in translation as a function of translation (rho/rho) [motionmodel]*/
    double getsrr() const { return m_motionModel.srr; }
    void setsrr(double srr) { m_motionModel.srr = srr; }

    /**odometry error in translation as a function of rotation (rho/theta) [motionmodel]*/
    double getsrt() const { return m_motionModel.srt; }
    void setsrt(double srt) { m_motionModel.srt = srt; }

    /**odometry error in rotation as a function of translation (theta/rho) [motionmodel]*/
    double getstr() const { return m_motionModel.str; }
    void setstr(double str) { m_motionModel.str = str; }

    /**odometry error in  rotation as a function of rotation (theta/theta) [motionmodel]*/
    double getstt() const { return m_motionModel.stt; }
    void setstt(double stt) { m_motionModel.stt = stt; }
		
    /**minimum score for considering the outcome of the scanmatching good*/
    double getminimumScore() const { return m_minimumScore; }
    void setminimumScore(double minimumScore) { m_minimumScore = minimumScore;}

  protected:
    /**Copy constructor*/
    GridSlamProcessor(const GridSlamProcessor& gsp);
 
    double m_minimumScore;
    double m_resampleThreshold;

    /**the laser beams*/
    unsigned int m_beams;
    double last_update_time_;
    double period_;
    
    
    /**the particles*/
    ParticleVector m_particles;

    /**the particle indexes after resampling (internally used)*/
    std::vector<unsigned int> m_indexes;

    /**the particle weights (internally used)*/
    std::vector<double> m_weights;
    
    /**the motion model*/
    MotionModel m_motionModel;

    /**this sets the neff based resampling threshold*/
  public:
    double getresampleThreshold() const { return m_resampleThreshold; }
    void setresampleThreshold(double resampleThreshold) { m_resampleThreshold = resampleThreshold; }
    double getneff() const { return m_neff; }

    //processing parameters (size of the map)
    double getxmin() const { return m_xmin; }
    double getymin() const { return m_ymin; }
    double getxmax() const { return m_xmax; }
    double getymax() const { return m_ymax; }

    //processing parameters (resolution of the map)
    double getdelta() const { return m_delta; }

    //registration score (if a scan score is above this threshold it is registered in the map)
    double getregScore() const { return m_regScore; }
    void setregScore(double regScore) { m_regScore = regScore; }

    //registration score (if a scan score is below this threshold a scan matching failure is reported)
    double getcritScore() const { return m_critScore; }
    void setcritScore(double critScore) { m_critScore = critScore; }    

    //registration score maximum move allowed between consecutive scans
    double getmaxMove() const { return m_maxMove; }
    void setmaxMove(double maxMove) { m_maxMove = maxMove; }

    //process a scan each time the robot translates of linearThresholdDistance
    double getlinearThresholdDistance() const { return m_linearThresholdDistance; }
    void setlinearThresholdDistance(double linearThresholdDistance) { m_linearThresholdDistance = linearThresholdDistance; }

    //process a scan each time the robot rotates more than angularThresholdDistance
    double getangularThresholdDistance() const { return m_angularThresholdDistance; }
    void setangularThresholdDistance(double angularThresholdDistance) { m_angularThresholdDistance = angularThresholdDistance; }

    //smoothing factor for the likelihood
    double getobsSigmaGain() const { return m_obsSigmaGain; }
    void setobsSigmaGain(double obsSigmaGain) { m_obsSigmaGain = obsSigmaGain; }

  protected:    
    //state
    int  m_count, m_readingCount;
    OrientedPoint m_lastPartPose;
    OrientedPoint m_odoPose;
    OrientedPoint m_pose;
    double m_linearDistance, m_angularDistance;
    double m_neff, m_delta, m_regScore, m_critScore;
    double m_xmin, m_ymin, m_xmax, m_ymax;
    double m_maxMove, m_linearThresholdDistance, m_angularThresholdDistance, m_obsSigmaGain;
	
    //stream in which to write the gfs file
    std::ofstream m_outputStream;

    // stream in which to write the messages
    std::ostream& m_infoStream;   
    
    // the functions below performs side effect on the internal structure,
    //should be called only inside the processScan method
  private:
    
    /**scanmatches all the particles*/
    inline void scanMatch(const double *plainReading);
    /**normalizes the particle weights*/
    inline void normalize();
    
    // return if a resampling occured or not
    inline bool resample(const double* plainReading, int adaptParticles, 
			 const RangeReading* rr=0);
    
    //tree utilities
    void updateTreeWeights(bool weightsAlreadyNormalized = false);
    void resetTree();
    double propagateWeights();
    
    // Messageing Functions.
    void printOdomJumpWarnings(const OrientedPoint & new_pose);

    void printBasicInfoForScanUpdate();

    void printUpdatedNeff();

    void printUpdatedParticlePose(const RangeReading & reading);

    void printScanMatchUpdateInfo(const RangeReading & reading);

    void printScanFrameInfo(const RangeReading & reading);

    void printFlush();

  };

typedef std::multimap<const GridSlamProcessor::TNode*, GridSlamProcessor::TNode*> TNodeMultimap;


#include "gridslamprocessor.hxx"

};

#endif
