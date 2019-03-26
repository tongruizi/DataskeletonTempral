#ifndef CORRECTSTRAIGHTENING_H
#define CORRECTSTRAIGHTENING_H

#include "StraighteningMethods.h"
#include "SuperSort.h"
#include "SegmentCompare.h"
#include "GeneralConvertor.h"


class CorrectStraightening
{
public:
    CorrectStraightening() {}
    //! First we do the allocation
//
//    static void GraphSplittingTools(MyGraphType & G)
//    {
//        std::deque<std::pair<vertex_descriptor,vertex_descriptor>> ddd;
//        auto vpair = boost::vertices(G);
//        for (auto it = vpair.first; it != vpair.second ; it++)
//        {
//            if (boost::degree(*it, G) != 2)
//            {
//                ddd.push_back(std::make_pair(*it,*it));
//                break;
//            }
//        }
//        while (!ddd.empty())
//        {
//            vertex_descriptor a = ddd.front().first;
//            vertex_descriptor wrong = ddd.front().second;
//
//            //     std::cout << a << std::endl;
//            ddd.pop_front();
//            auto epair = boost::adjacent_vertices(a,G);
//            for (auto it = epair.first ; it != epair.second; it++)
//            {
//                //  std::cout << *it << std::endl;
//                if (*it == wrong)
//                {
//                    continue;
//                }
//                std::vector<PointPlus> path;
//                addToPath(path,G,a);
//                vertex_descriptor ot = *it;
//                vertex_descriptor prev = a;
//                while (boost::degree(ot,G) == 2)
//                {
//                    //    std::cout << "InnerLoop: " << ot << std::endl;
//                    auto innerpair = boost::adjacent_vertices(ot,G);
//                    for (auto bt = innerpair.first ; bt != innerpair.second; bt++)
//                    {
//                        if (*bt != prev)
//                        {
//                            prev = ot;
//                            ot = *bt;
//                            break;
//                        }
//                    }
//                    addToPath(path,G,ot);
//
//                }
//                addToPath(path,G,ot);
//                paths.push_back(path);
//                if (boost::degree(ot,G) != 1)
//                {
//                    ddd.push_back(std::make_pair(ot, prev));
//                }
//            }
//        }
//
//    }

    static void InitialAllocation(MyGraphType & G, std::vector<int> & switchIndices, std::vector<Segment> & theSegments)
    {
        std::vector<std::vector<PointPlus>> temppaths;
        StraighteningMethods::GraphToPaths(temppaths,G);
        //! Convert temppaths to theSegments
        for (int i = 0; i < temppaths.size(); i++)
        {
            switchIndices.push_back(theSegments.size());
            for (int j = 0; j < temppaths[i].size() - 1; j++)
            {
                theSegments.push_back(Segment(temppaths[i][j].first,temppaths[i][j+1].first));
            }
        }
    //    std::cout << "Initial allocation suceful" << std::endl;

    }


    static void AllocationOfSegments(std::vector<int> & switchIndices, std::vector<Segment> & theSegments, std::list<Point> & cloud,
                                     std::vector<std::vector<Point>> & theAllocatedPoints)
    {
      //  std::cout << "Allocating segments..." << std::endl;
        std::map<Segment,int,SegmentCompare> theMap;
        for (int i = 0; i < theSegments.size(); i++)
        {
            theMap[theSegments[i]] = i;
        }
     //   std::cout << "Stil lfine..." << std::endl;

        VectorSegmentTree AABB(theSegments.begin(),theSegments.end());
        AABB.accelerate_distance_queries();
        theAllocatedPoints.resize(theSegments.size());

      //  std::cout << "AABB tree good..." << std::endl;

        for (auto it = cloud.begin(); it != cloud.end(); it++)
        {
            auto qwe = AABB.closest_point_and_primitive(*it);
            // std::cout << "The index:" << theMap[*(qwe.second)] << std::endl;
            theAllocatedPoints[theMap[*(qwe.second)]].push_back(*it);
        }
        for (int i = 0; i < theAllocatedPoints.size(); i++)
        {
            Point p = theSegments[i].source();
            std::sort(theAllocatedPoints[i].begin(), theAllocatedPoints[i].end(),[p](Point & a, Point & b)
            {
                return (CGAL::squared_distance(p,a)) < (CGAL::squared_distance(p,b));
            });
            theAllocatedPoints[i].insert(theAllocatedPoints[i].begin(),theSegments[i].source());
            theAllocatedPoints[i].push_back(theSegments[i].target());

        }
      //  std::cout << "Exited allocation of segments method" << std::endl;
    }

    static double maximalDistance(int i, int j, std::vector<Point> & points)
    {
        Segment s = Segment(points[i],points[j]);
        double maxD = 0;
        for (int k = i; k <= j; k++)
        {
            maxD = std::max(maxD, sqrt(CGAL::squared_distance(points[k],s)));
        }
    }


    static bool RunTheAlgorithm(std::vector<Point> & points, std::list<Point> & thepath, double e)
    {
        int sizepath = points.size();
        int a = 0;
        int b = 0;
        while (b < points.size() - 1)
        {
            int l = 0;
          //  l = l + 1; //! Is this required? Most likely now
            while((CorrectStraightening::maximalDistance(b, b + pow(2,(l+1)), points) <= e)
                    &&(b + pow(2,(l+1)) < points.size()))
            {
                l++;
            }
            int low = pow(2,l);
            int sizesizesize = points.size();
            int powerpowerpower = pow(2,l+1);
            int high = std::min(powerpowerpower,sizesizesize - b);
            while(low < high - 1)
            {
                int mid = (low + high)/2;
                int sumdxd = b + mid;
                double maxd = CorrectStraightening::maximalDistance(b, b + mid, points);
                if (maxd <= e)
                {
                    low = mid;
                }
                else
                {
                    high = mid;
                }
            }
            b = b + low;
            if (b < points.size())
            {
                thepath.push_back(points[b]);
            }
            else
            {
                thepath.push_back(points[points.size()-1]);
            }
            a = a + 1;
        }
        return true;
    }

    static void LaunchMechanism(int start, int endd,  std::vector<std::vector<Point>> & theAllocatedPoints, std::list<std::list<Point>> & out,
                                std::vector<Segment> & sss, double ce, std::vector<std::vector<Point>> & debuglist)
    {
      // std::cout << "Enetering launchmechanism: "  << " Start: " << start << " End: " << endd << std::endl;
      //  std::cout << "THEAllocatedPoints size: " << theAllocatedPoints.size() << std::endl;
        std::vector<Point> pointsConsidered;
        std::list<Point> finalpath;
        Point firstPoint = sss[start].source();
        Point lastPoint = sss[endd-1].target();
        finalpath.push_back(firstPoint);
        for (int i = start ; i < endd; i++)
        {
            for (int q = 0; q < theAllocatedPoints[i].size(); q++)
            {
                pointsConsidered.push_back(theAllocatedPoints[i][q]);
            }

        }
      //  debuglist.push_back(pointsConsidered);
       // std::cout << "Size of points considered: " << pointsConsidered.size() << std::endl;
        pointsConsidered[0] = firstPoint;
        pointsConsidered[pointsConsidered.size()-1] = lastPoint;
        //! Run the algorithm:
      // std::cout << "Before running the algorithm: " << " Start: " << start << " End: " << endd << std::endl;
        RunTheAlgorithm(pointsConsidered,finalpath,ce);
     //  std::cout << "the main algorithm finnished runnning " << std::endl;
        if (finalpath.front() == finalpath.back())
        {
        std::cout << "Adding bad boy" << std::endl;
        std::cout << "!!!!!!!!!!!!!!!" << std::endl;
         std::cout << "!!!!!!!!!!!!!!!" << std::endl;
          std::cout << "!!!!!!!!!!!!!!!" << std::endl;
        }
        out.push_back(finalpath);
    }

    static double ComputeStraightening(MyGraphType & G, std::list<Point> & cloud, std::list<std::list<Point>> & out, double e)
    {
    //    std::cout << "Method accessed" << std::endl;
        std::vector<std::vector<Point>> debuglist;
        std::vector<int> switchIndices;
        std::vector<Segment> theSegments;
        std::vector<std::vector<Point>> theAllocatedPoints;
        CorrectStraightening::InitialAllocation(G,  switchIndices, theSegments);
        CorrectStraightening::AllocationOfSegments(switchIndices, theSegments,  cloud, theAllocatedPoints);
     //   std::cout << "Preparing to compute ce: " << std::endl;
   //     std::cout << "Size of segments: " << theSegments.size() << std::endl;
    //    std::cout << "TheAllocatedPoints size: " << theAllocatedPoints.size() << std::endl;
        int sumsize = 0;
        for (int i = 0; i < theAllocatedPoints.size(); i++)
        {
        //    std::cout << "Elements in: " << theAllocatedPoints[i].size() << std::endl;
            sumsize = sumsize + theAllocatedPoints[i].size();
        }
        double ce = 0;
        //! Coompute ce:
        for (int i = 0; i < theSegments.size(); i++)
        {
            for (int j = 0; j < theAllocatedPoints[i].size(); j++)
            {
                ce = std::max(ce, sqrt(CGAL::squared_distance(theSegments[i],theAllocatedPoints[i][j])));

            }
        }
        double cr = ce * e;
     //   std::cout << "Initial error is computed to be: " << ce << std::endl;
 //      std::cout << "Stuff allocated" << std::endl;
        //! Prepare the correctStraightening:

        int start = 0;
        int endd = 0;
     //   std::cout << "SWitchIndices size: " << switchIndices.size() << std::endl;
        for (int www = 0; www < switchIndices.size() - 1; www++ )
        {
            start = switchIndices[www];
            endd = switchIndices[www+1];
            LaunchMechanism(start,endd,theAllocatedPoints,out,theSegments,cr,debuglist);
        }
      //  std::cout << "Exited the first loop" << std::endl;
        LaunchMechanism(switchIndices[switchIndices.size()-1],theSegments.size(), theAllocatedPoints,out,theSegments,cr,debuglist);
     //   std::cout << "Method exited succefully" << std::endl;
    //   std::cout << "Number of paths: " << out.size() << std::endl;
     //   for  (auto it = out.begin(); it != out.end(); it++)
     //   {
    //        std::cout << (*it).size() << std::endl;

    //    }
   //     GeneralConvertor::StraighteningDebugPrint("/home/yury/LocalTests/RealCloudDebug2/debug.csv",debuglist);
   //     std::cout << "Value is: " << ce << std::endl;
   //     std::cout << "Sumsize: " << sumsize << std::endl;
        return ce;
    }

protected:

private:
};

#endif // CORRECTSTRAIGHTENING_H
