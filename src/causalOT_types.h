#include <RcppEigen.h>

typedef Eigen::VectorXd vector;
typedef Eigen::Matrix<long double, Eigen::Dynamic,  1> vectorLD;
typedef Eigen::VectorXi vectorI;
typedef Eigen::RowVectorXd rowVector;

typedef Eigen::MatrixXd matrix;
typedef Eigen::MatrixXi matrixI;
typedef Eigen::Matrix<long double, Eigen::Dynamic,  Eigen::Dynamic> matrixLD;

typedef Eigen::Ref<matrix> refMat;
typedef Eigen::Ref<matrixI> refMatI;

typedef Eigen::Ref<Eigen::ArrayXi> refArrayI;
typedef Eigen::Ref<Eigen::ArrayXd> refArray;

typedef Eigen::Ref<vector> refVec;
typedef Eigen::Ref<vectorI> refVecI;

typedef Eigen::Ref<const matrix> refMatConst;
typedef Eigen::Ref<const matrixI> refMatConstI;
typedef Eigen::Ref<const Eigen::ArrayXd> refArrayConst;
typedef Eigen::Ref<const Eigen::ArrayXi> refArrayConstI;

typedef Eigen::Ref<const vector> refVecConst;
typedef Eigen::Ref<const vectorI> refVecConstI;

typedef matrix::ColXpr ColXpr;
typedef matrixI::ColXpr ColXprI;

typedef matrix::RowXpr RowXpr;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rowMat;
typedef Eigen::LLT<matrix> llt;
typedef Eigen::LDLT<matrix> ldlt;

typedef Eigen::Map<matrix> matMap;
typedef Eigen::Map<const matrix> matMapConst;
typedef Eigen::Map<matrixLD> matMapLD;
typedef Eigen::Map<rowMat> rowMatMap;

typedef Eigen::Map<Eigen::VectorXd> vecMap;
typedef Eigen::Map<vectorLD> vecMapLD;
typedef Eigen::Map<const vector> vecMapConst;
typedef Eigen::Map<Eigen::VectorXi> vecMapI;
typedef Eigen::Map<const Eigen::VectorXi> vecMapConstI;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseVector<double> SpVec;

