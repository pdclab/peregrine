#ifndef OPTIONS_HH
#define OPTIONS_HH

namespace Peregrine
{

  /**
   * Indicates when aggregation occurs:
   * either AT_THE_END of execution or ON_THE_FLY
   */
  enum OnTheFlyOption
  {
    ON_THE_FLY,
    AT_THE_END,
  };
  
  /**
   * Toggle for early termination.
   */
  enum StoppableOption
  {
    STOPPABLE,
    UNSTOPPABLE,
  };
}

#endif
