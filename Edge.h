#ifndef EDGE_H
#define EDGE_H

#include "Node.h"

class Edge{
    public:
          //Edge();
          Edge(std::string edgeName, Node a, Node b); // constructor
          ~Edge();
          Node getHead()const;
          Node getTail()const;
          Node getmidpoint()const;
          float getLength()const;
          std::string getEdgeName()const;
          float getXComponentNorm()const;
          float getYComponentNorm ()const;
        //  Node *getNorme()const;
          //float * getNorme()const;
          float getradius()const;
          bool operator==(const Edge &edge) const {
          return edge.name== name;
          }


    private:
            Node head;
            Node tail;
            float length;
            //Node *midpoint;
            //Node *norme;
            std::string name;
            float XComponentNorm;
            float YComponentNorm;
            float radius;
            friend std::ostream& operator<<(std::ostream&, const Edge&);
};
#endif
