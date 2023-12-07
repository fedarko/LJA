#pragma once

#include "assembly_graph.hpp"
#include "sequences/contigs.hpp"

namespace ag {
    template<class Traits>
    class GraphPath {
    public:
        typedef typename Traits::Vertex Vertex;
        typedef typename Traits::Edge Edge;
    private:
        Vertex *start_;
        std::vector<Edge *> path;
        size_t cut_left = 0;
        size_t cut_right = 0;
    public:
        typedef typename std::vector<Edge *>::iterator iterator;
        typedef typename std::vector<Edge *>::const_iterator const_iterator;
        typedef TransformingIterator<CountingIterator<size_t>, Vertex> vertex_iterator;
        typedef TransformingIterator<CountingIterator<size_t>, Edge> edge_iterator;
        typedef TransformingGenerator<CountingIterator<size_t>, Segment<Edge>> segment_iterator;

        GraphPath(Vertex &_start, std::vector<Edge *> _path, size_t left_skip, size_t rightSkip) :
                    start_(&_start), path(std::move(_path)), cut_left(left_skip), cut_right(rightSkip) {}
        GraphPath(Vertex &_start) : start_(&_start) {} // NOLINT(google-explicit-constructor)
        GraphPath(Vertex &_start, size_t left_skip, size_t right_skip) : start_(&_start), cut_left(left_skip), cut_right(right_skip) {} // NOLINT(google-explicit-constructor)
        GraphPath(Edge &edge) : start_(&edge.getStart()), path({&edge}) {} // NOLINT(google-explicit-constructor)
        GraphPath(const Segment<Edge> &segment) : start_(&segment.contig().getStart()), // NOLINT(google-explicit-constructor)
                          path({&segment.contig()}), cut_left(segment.left), cut_right(segment.contig().truncSize() - segment.right) {}
        GraphPath() : start_(nullptr) {}
        template<class Iterator>
        explicit GraphPath(Iterator begin, Iterator end) : cut_left(0), cut_right(0) {
            while(begin != end) {
                *this += *begin;
                ++begin;
            }
            if(size() > 0) {
                start_ = &frontEdge().getStart();
            }
        }

        static GraphPath WalkForward(Edge &start);

        Vertex &getVertex(size_t i) const;
        Edge &getEdge(size_t i) const {return *path[i];};
        Vertex &start() const {return *start_;}
        Vertex &finish() const {return path.empty() ? start() : path.back()->getFinish();}
        size_t find(Edge &edge, size_t pos = 0) const;
        size_t find(Vertex &v, size_t pos = 0) const;
        Edge &backEdge() const {return *path.back();}
        Edge &frontEdge() const {return *path.front();}
        size_t size() const {return path.size();}
        bool empty() const {return path.empty();}

    //        TODO: Find a way to iterate over temporary path objects
        IterableStorage<vertex_iterator> vertices() const &;
        IterableStorage<vertex_iterator> vertices() && = delete;
        IterableStorage<edge_iterator> edges() const &;
        IterableStorage<edge_iterator> edges() && = delete;

        segment_iterator begin() const;
        segment_iterator  end() const;

        GraphPath RC() const;
        double minCoverage() const;
        Sequence Seq() const;
        Sequence truncSeq() const;
        size_t truncLen() const;
        size_t len() const;

        Segment<Edge> back() const;
        Segment<Edge> front() const;
        Segment<Edge> operator[](size_t i) const;
        std::string covStr(bool show_coverage = false) const;
        std::string str() const;
        std::string lenStr() const;

        bool valid() const;
        void invalidate();
        GraphPath subPath(size_t from, size_t to) const;
        GraphPath subPath(size_t from) const {return subPath(from, size());}
        GraphPath reroute(size_t left, size_t right, const GraphPath &rerouting) const;
        void operator+=(const GraphPath &other);
        void operator+=(const Segment<Edge> &other);
        void operator+=(Edge &other);
        GraphPath operator+(const GraphPath &other) const;
        GraphPath operator+(const Segment<Edge> &other) const;
        GraphPath operator+(Edge &other) const;
        //TODO deprecate
        Sequence map(std::unordered_map<const Edge *, Sequence> &edge_map);

        void pop_back() {
            path.pop_back();
            cut_right = 0;
            if(path.empty() && cut_left > 0) {
                invalidate();
            }
        }
        void pop_back(size_t len) {
            path.erase(path.end() - len, path.end());
            if(len > 0)
                cut_right = 0;
            if(size() == 0 && cut_left != 0)
                invalidate();
        }
        GraphPath &cutBack(size_t l);
        GraphPath &cutFront(size_t l);
        GraphPath &uniqueExtendBack(size_t l);
        GraphPath &uniqueExtendFront(size_t l);
        GraphPath &addStep();
        GraphPath &addStep(Edge &edge);
        std::vector<GraphPath> allSteps();
        std::vector<GraphPath> allExtensions(size_t len);
        GraphPath &extend(const Sequence &seq);

        unsigned char lastNucl() const;
        size_t leftSkip() const;
        size_t rightSkip() const;
        bool endClosed() const;
        bool startClosed() const;

        Sequence truncSubseq(size_t start_position, size_t size = 10000000) const;

        bool operator==(const GraphPath &other) const {
            return start_ == other.start_ && cut_left == other.cut_left && cut_right == other.cut_right && path == other.path;
        }
        bool operator!=(const GraphPath &other) const {return !operator==(other);}
    };
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::WalkForward(Edge &start) {
    GraphPath<Graph> res(start);
    Vertex *next = &start.getFinish();
    VERIFY(next != nullptr);
    while(*next != start.getStart() && *next != start.getStart().rc() && !next->isJunction()) {
        VERIFY(next != nullptr);
        res += next->front();
        next = &res.finish();
    }
    return std::move(res);
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::subPath(size_t from, size_t to) const {
    if (from == to) {
        if ((from == 0 && leftSkip() > 0) || (to == size() && rightSkip() > 0))
            return {};
        else
            return GraphPath<Graph>(getVertex(from));
    } else
        return {getVertex(from), std::vector<Edge *>(path.begin() + from, path.begin() + to),
                from == 0 ? leftSkip() : 0, to == size() ? rightSkip() : 0};
}

template<class Graph>
typename Graph::Vertex &ag::GraphPath<Graph>::getVertex(size_t i) const {
    VERIFY(i <= path.size());
    if (i == 0)
        return *start_;
    else
        return path[i - 1]->getFinish();
}

template<class Graph>
size_t ag::GraphPath<Graph>::find(Edge &edge, size_t pos) const {
    while(pos < size() && edge != *path[pos])
        pos++;
    if(pos == size())
        return -1;
    return pos;
}

template<class Graph>
size_t ag::GraphPath<Graph>::find(Vertex &v, size_t pos) const {
    while(pos <= size() && v != getVertex(pos))
        pos++;
    if(pos > size())
        return -1;
    return pos;
}

template<class Graph>
double ag::GraphPath<Graph>::minCoverage() const {
    double res = 100000000;
    for (const Edge *edge : path) {
        res = std::min(edge->getCoverage(), res);
    }
    return res;
}

template<class Graph>
size_t ag::GraphPath<Graph>::truncLen() const {
    size_t res = 0;
    for (Edge *edge : path)
        res += edge->truncSize();
    if(res == 0)
        return 0;
    return res - cut_left - cut_right;
}

template<class Graph>
IterableStorage<typename ag::GraphPath<Graph>::vertex_iterator> ag::GraphPath<Graph>::vertices() const & {
    std::function<Vertex &(size_t)> transformer = [this](size_t ind)->Vertex &{return getVertex(ind);};
    CountingIterator<size_t> end_it = CountingIterator<size_t>(valid()? size() + 1 : 0);
    vertex_iterator vbegin (CountingIterator<size_t>(0), end_it, transformer);
    vertex_iterator vend(end_it, end_it, transformer);
    return {vbegin, vend};
}

template<class Graph>
IterableStorage<typename ag::GraphPath<Graph>::edge_iterator> ag::GraphPath<Graph>::edges() const &{
    std::function<Edge &(size_t)> transformer = [this](size_t ind)->Edge &{return *path[ind];};
    CountingIterator<size_t> end_it = CountingIterator<size_t>(size());
    edge_iterator ebegin (CountingIterator<size_t>(0), end_it, transformer);
    edge_iterator eend(end_it, end_it, transformer);
    return {ebegin, eend};
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::RC() const {
    if(!valid())
        return {};
    std::vector<Edge *> res;
    for(auto it  = path.rbegin(); it != path.rend(); ++it) {
        Edge *e = *it;
        res.emplace_back(&e->rc());
    }
    return {finish().rc(), res, rightSkip(), leftSkip()};
}

template<class Graph>
void ag::GraphPath<Graph>::invalidate() {
    start_ = nullptr;
    path.clear();
    cut_left = 0;
    cut_right = 0;
}

template<class Graph>
bool ag::GraphPath<Graph>::valid() const {
    VERIFY(start_ != nullptr || size() == 0);
    return start_ != nullptr;
}

template<class Graph>
ag::GraphPath<Graph> &ag::GraphPath<Graph>::cutBack(size_t l) {
    VERIFY(l <= len());
    size_t expected = len() - l;
    size_t cur_cut = 0;
    size_t cut = 0;
    l += cut_right;
    cut_right = 0;
    while(cur_cut < size() && l >= path[path.size() - 1 - cur_cut]->truncSize()) {
        if(path[path.size() - 1 - cur_cut]->truncSize() == 0) {
            cur_cut++;
        } else {
            l -= path[path.size() - 1 - cur_cut]->truncSize();
            cur_cut++;
            cut = cur_cut;
        }
    }
    if(cut == size()) {
        *this = {start(), leftSkip(), l};
    } else {
        path.erase(path.end() - cut, path.end());
        cut_right = l;
    }
    VERIFY(len() == expected);
    return *this;
}

template<class Graph>
ag::GraphPath<Graph> &ag::GraphPath<Graph>::uniqueExtendBack(size_t l) {
    if(cut_right != 0) {
        size_t tmp = std::min(l, cut_right);
        l -= tmp;
        cut_right -= tmp;
    }
    while(l > 0) {
        VERIFY(finish().outDeg() == 1);
        Edge &e = finish().front();
        size_t tmp = std::min(e.truncSize(), l);
        *this += e;
        cutBack(e.truncSize() - tmp);
    }
    return *this;
}


//TODO: Optimize
template<class Graph>
ag::GraphPath<Graph> &ag::GraphPath<Graph>::uniqueExtendFront(size_t l) {
    *this = this->RC().uniqueExtendBack(l);
    return *this;
}

template<class Graph>
ag::GraphPath<Graph> &ag::GraphPath<Graph>::cutFront(size_t l) {
    VERIFY(l <= len());
    size_t expected = len() - l;
    size_t cur_cut = 0;
    size_t cut = 0;
    l += cut_left;
    cut_left = 0;
    while(cur_cut < size() && l >= path[cur_cut]->rc().truncSize()) {
        if(path[cur_cut]->rc().truncSize() == 0) {
            cur_cut++;
        } else {
            l -= path[cur_cut]->rc().truncSize();
            cur_cut++;
            cut = cur_cut;
        }
    }
    if(cut == size()) {
        *this = {finish(), l, rightSkip()};
    } else {
        path.erase(path.begin(), path.begin() + cut);
        start_ = &path.front()->getStart();
        cut_left = l;
    }
    VERIFY(len() == expected);
    return *this;
}

template<class Graph>
ag::GraphPath<Graph> &ag::GraphPath<Graph>::addStep() {
    cut_right -= 1;
    return *this;
}

template<class Graph>
ag::GraphPath<Graph> &ag::GraphPath<Graph>::addStep(Edge &edge) {
    *this += Segment<Edge>(edge, 0, 1);
    return *this;
}

template<class Graph>
ag::GraphPath<Graph> &ag::GraphPath<Graph>::extend(const Sequence &seq) {
    VERIFY(valid());
    for (size_t cpos = 0; cpos < seq.size(); cpos++) {
        unsigned char c = seq[cpos];
        if (endClosed()) {
            Vertex &v = finish();
            if (v.hasOutgoing(c)) {
                Edge &edge = v.getOutgoing(c);
                addStep(edge);
            } else {
                invalidate();
                return *this;
            }
        } else {
            if (back().contig().truncSeq()[back().right] == c) {
                addStep();
            } else {
                invalidate();
                return *this;
            }
        }
    }
    return *this;
}

template<class Graph>
bool ag::GraphPath<Graph>::endClosed() const {
    return valid() && rightSkip() == 0;
}

template<class Graph>
bool ag::GraphPath<Graph>::startClosed() const {
    return valid() && leftSkip() == 0;
}

template<class Graph>
unsigned char ag::GraphPath<Graph>::lastNucl() const {
    Segment<Edge> seg = back();
    return seg.truncSeq()[seg.right - 1];
}

template<class Graph>
size_t ag::GraphPath<Graph>::leftSkip() const {
    return cut_left;
}

template<class Graph>
size_t ag::GraphPath<Graph>::rightSkip() const {
    return cut_right;
}

template<class Graph>
std::vector<ag::GraphPath<Graph>> ag::GraphPath<Graph>::allSteps() {
    if (size() != 0 && cut_right > 0) {
        ag::GraphPath<Graph> copy = *this;
        return {std::move(copy.addStep())};
    }
    std::vector<ag::GraphPath<Graph>> res;
    Vertex &end = finish();
    for (Edge &edge : end) {
        GraphPath<Graph> copy = *this;
        res.emplace_back(std::move(copy.addStep(edge)));
    }
    return res;
}

template<class Graph>
std::vector<ag::GraphPath<Graph>> ag::GraphPath<Graph>::allExtensions(size_t len) {
    std::vector<GraphPath<Graph>> res = {*this};
    size_t left = 0;
    size_t right = 1;
    for (size_t l = 0; l < len; l++) {
        for (size_t i = left; i < right; i++) {
            std::vector<GraphPath<Graph>> tmp = res[i].allSteps();
            res.insert(res.end(), tmp.begin(), tmp.end());
        }
        left = right;
        right = res.size();
    }
    return std::move(res);
}

template<class Graph>
Sequence ag::GraphPath<Graph>::map(std::unordered_map<const Edge *, Sequence> &edge_map) {
    SequenceBuilder sb;
    bool start = true;
    for (Segment<Edge> seg : *this) {
        auto it = edge_map.find(&seg.contig());
        if (it == edge_map.end()) {
            if (start) {
                sb.append((start_->getSeq() + seg.contig().truncSeq()).Subseq(seg.left, seg.right + start_->getSeq().size()));
                start = false;
            } else {
                sb.append(seg.truncSeq());
            }
        } else {
            size_t left = start_->getSeq().size();
            if (start) {
                left = 0;
            }
            size_t right = start_->getSeq().size();
            size_t sz = it->second.size() - start_->getSeq().size();
            if (seg.left == 0 && seg.right == seg.contig().size()) {
                right += sz;
            } else if (seg.left == 0) {
                right += std::min(sz, seg.right);
            } else if (seg.right == seg.contig().size()) {
                left += sz - std::min(sz, seg.size());
                right += sz;
            } else {
                size_t l = seg.left * sz / seg.contig().size();
                left += l;
                right += std::min(l + seg.size(), sz);
            }
            sb.append(it->second.Subseq(left, right));
            start = false;
        }
    }
    return sb.BuildSequence();
}

template<class Graph>
Sequence ag::GraphPath<Graph>::Seq() const {
    if (!valid())
        return {};
    if(size() == 0) {
        return start_->getSeq().Subseq(leftSkip(), start_->size() - rightSkip());
    }
    Edge & edge = frontEdge();
    Sequence seq = edge.getSeq();
    if(size() == 1) {
        return seq.Subseq(leftSkip(), seq .size() - rightSkip());
    }

    SequenceBuilder sb;
    sb.append(seq.Subseq(leftSkip()));
    for(size_t i = 1; i + 1 < size(); i++) {
        sb.append(operator[](i).truncSeq());
    }
    sb.append(backEdge().truncSeq().Subseq(0, backEdge().truncSize() - rightSkip()));
    return sb.BuildSequence();
}

template<class Graph>
Sequence ag::GraphPath<Graph>::truncSeq() const {
    SequenceBuilder sb;
    for (Segment<Edge> seg : *this) {
        sb.append(seg.truncSeq());
    }
    return sb.BuildSequence();
}

template<class Graph>
Sequence ag::GraphPath<Graph>::truncSubseq(size_t start_position, size_t sz) const {
    SequenceBuilder sb;
    for (size_t i = start_position; i < this->size(); i++) {
//            std::cout << i << " " << sz << " " << size << std::endl;
//            std::cout << als[i].contig().size() << " " << als[i].left << " " << als[i].right << " " << als[i].size() << std::endl;
        Segment<Edge> seg = (*this)[i];
        if (seg.size() >= sz) {
            sb.append(seg.shrinkRightToLen(sz).truncSeq());
            sz = 0;
            break;
        } else {
            sb.append(seg.truncSeq());
            sz -= seg.size();
        }
    }
    return sb.BuildSequence();
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::reroute(size_t left, size_t right, const GraphPath<Graph> &rerouting) const {
    VERIFY(left == 0 || getVertex(left) == rerouting.start());
    VERIFY(right == size() || getVertex(right) == rerouting.finish());
    GraphPath<Graph> res;
    res += subPath(0, left);
    res += rerouting;
    res += subPath(right, size());
    return std::move(res);
}

template<class Graph>
void ag::GraphPath<Graph>::operator+=(const ag::GraphPath<Graph> &other) {
    if(other.size() == 0)
        return;
    if(!valid()) {
        *this = other;
        return;
    }
    VERIFY(finish() == other.getVertex(0));
    for (Segment<Edge> al : other) {
        operator+=(al);
    }
}

template<class Graph>
void ag::GraphPath<Graph>::operator+=(const Segment<Edge> &other) {
    if(!valid()) {
        *this = {other};
        return;
    }
    if(cut_right == 0) {
        VERIFY(other.left == 0 && finish() == other.contig().getStart());
        path.emplace_back(&other.contig());
        cut_right = other.contig().truncSize() - other.right;
    } else {
        VERIFY(cut_right == other.contig().truncSize() - other.left && other.contig() == backEdge());
    }
    cut_right = other.contig().truncSize() - other.right;
}

template<class Graph>
void ag::GraphPath<Graph>::operator+=(Edge &other) {
    GraphPath<Graph>::operator+=(Segment<Edge>(other, 0, other.truncSize()));
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::operator+(const GraphPath<Graph> &other) const {
    GraphPath<Graph> res = *this;
    res += other;
    return std::move(res);
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::operator+(const Segment<Edge> &other) const {
    GraphPath<Graph> res = *this;
    res += other;
    return std::move(res);
}

template<class Graph>
ag::GraphPath<Graph> ag::GraphPath<Graph>::operator+(Edge &other) const {
    GraphPath<Graph> res = *this;
    res += other;
    return std::move(res);
}

template<class Graph>
std::string ag::GraphPath<Graph>::covStr(bool show_coverage) const {
    if(!valid())
        return "";
    std::stringstream ss;
    ss << leftSkip() << " " << start().getInnerId();
    for(const Segment<Edge> &seg : *this) {
        ss << " " << seg.size() << "/" <<seg.contig().truncSize() << seg.contig().nuclLabel() ;
        if(show_coverage) {
            ss << "(" << seg.contig().getCoverage() << ")";
        }
        ss << " " << seg.contig().getFinish().getInnerId();
    }
    ss << " " << rightSkip();
    return ss.str();
}

template<class Graph>
std::string ag::GraphPath<Graph>::str() const {
    if(!valid())
        return "";
    std::stringstream ss;
    ss << leftSkip() << " " << start().getInnerId();
    for(const Segment<Edge> &seg : *this) {
        ss << " " << seg.size() << "/" <<seg.contig().truncSize() << seg.contig().nuclLabel() ;
        ss << " " << seg.contig().getFinish().getInnerId();
    }
    ss << " " << rightSkip();
    return ss.str();
}

template<class Graph>
std::string ag::GraphPath<Graph>::lenStr() const {
    if(!valid())
        return "";
    std::stringstream ss;
    ss << leftSkip() << " [" << start().getInnerId() << "(" << start().size() << ")";
    for(const typename Graph::Edge &edge : edges()) {
        ss << " -> " << edge.nuclLabel() << "(" << edge.fullSize() << ") -> "  << edge.getFinish().getInnerId() << "(" << edge.getFinish().size() << ")" ;
    }
    ss << "] " << rightSkip();
    return ss.str();
}

template<class Graph>
Segment<typename Graph::Edge> ag::GraphPath<Graph>::back() const {
    return {*path.back(), (size() == 1 ? leftSkip() : 0), path.back()->truncSize() - rightSkip()};
}

template<class Graph>
Segment<typename Graph::Edge> ag::GraphPath<Graph>::front() const {
    return {*path.front(), leftSkip(), size() == 1 ? path.front()->truncSize() - rightSkip() : path.front()->truncSize()};
}

template<class Graph>
Segment<typename Graph::Edge> ag::GraphPath<Graph>::operator[](size_t i) const {
    return {*path[i], i == 0 ? leftSkip() : 0, i == size() - 1 ? path.back()->truncSize() - rightSkip() : path[i]->truncSize()};
}

template<class Graph>
typename ag::GraphPath<Graph>::segment_iterator ag::GraphPath<Graph>::begin() const {
    std::function<Segment<Edge> (size_t)> transformer = [this](size_t ind)->Segment<Edge> {
        return operator[](ind);
    };
    CountingIterator<size_t> end_it(size());
    return {CountingIterator<size_t>(0), end_it, transformer};
}

template<class Graph>
typename ag::GraphPath<Graph>::segment_iterator ag::GraphPath<Graph>::end() const {
    std::function<Segment<Edge> (size_t)> transformer = [this](size_t ind)->Segment<Edge> {
        return operator[](ind);
    };
    CountingIterator<size_t> end_it(size());
    return {end_it, end_it, transformer};
}

template<class Graph>
size_t ag::GraphPath<Graph>::len() const {
    size_t res = start().size();
    for(Edge &edge : edges()) {
        res += edge.truncSize();
    }
    return res - leftSkip() - rightSkip();
}
