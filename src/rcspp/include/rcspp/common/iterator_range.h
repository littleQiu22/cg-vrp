namespace rcspp{

template<typename Iterator>
class IteratorRange{
public:
    IteratorRange(Iterator begin, Iterator end): begin_(begin), end_(end){}
    Iterator begin() const {return begin_;}
    Iterator end() const {return end_;} 

private:
    Iterator begin_, end_;
};

}