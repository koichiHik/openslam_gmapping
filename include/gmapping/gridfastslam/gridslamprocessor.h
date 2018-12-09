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

    struct GridSlamProcParams {
      GridSlamProcParams()
        : particle_num(0), 
          period(5.0),
          obsSigmaGain(1.0),
          resampleThreshold(0.5),
          minimumScore(0.0)
      {}

      unsigned int particle_num;
      double period;
      double x_min, y_min, x_max, y_max;
      double delta;
      double maxMove;
      double linearThresholdDistance;
      double angularThresholdDistance;
      double resampleThreshold;
      double obsSigmaGain;
      double regScore, critScore;
      double minimumScore;
    };

    struct Params {
      GridSlamProcessor::GridSlamProcParams gridSlamProcParams;
      ScanMatcher::Params scanMatcherParams;
      MotionModel::Params motionModelParams;
    };
    
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
    
    /**This class defines a particle of the filter. Each particle has a map, a pose, 
     * a weight and retains the current node in the trajectory tree
     **/
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

    void init(const GridSlamProcessor::Params& params, OrientedPoint initialPose=OrientedPoint(0,0,0));

    void setGenerateMap(bool val) { m_matcher.setGenerateMap(val); }

    bool getGenerateMap() const { m_matcher.getGenerateMap(); }

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

  protected:
    /**Copy constructor*/
    GridSlamProcessor(const GridSlamProcessor& gsp);

    /**the laser beams*/
    unsigned int m_beams;
    double last_update_time_;  
    
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
    double getneff() const { return m_neff; }

    //processing parameters (size of the map)
    double getxmin() const { return m_xmin; }
    double getymin() const { return m_ymin; }
    double getxmax() const { return m_xmax; }
    double getymax() const { return m_ymax; }

    //processing parameters (resolution of the map)
    double getdelta() const { return m_delta; }

  protected:    
    //state
    int  m_count, m_readingCount;
    OrientedPoint m_lastPartPose;
    OrientedPoint m_odoPose;
    OrientedPoint m_pose;
    double m_neff;
    double m_xmin, m_ymin, m_xmax, m_ymax, m_delta;
    double m_linearDistance, m_angularDistance;
	
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

    Params m_params;

  };

typedef std::multimap<const GridSlamProcessor::TNode*, GridSlamProcessor::TNode*> TNodeMultimap;


#include "gridslamprocessor.hxx"

};

#endif
