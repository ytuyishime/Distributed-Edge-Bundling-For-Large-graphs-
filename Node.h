#ifndef NODE_H
#define NODE_H

//using namespace std;

class Node{
      public:
            Node();
            Node(std::string indexname, float x, float y);
            ~Node();
            float getX()const;
            float getY()const;
            void setX(float x);
            void setY(float y);
            bool issame(Node &p)const;
            void setIndex(std::string index);
            std::string getIndex()const;
      private:
        float x,y;
        std::string index;
        friend std::ostream& operator<<(std::ostream&, const Node&);
};
#endif
