// NodeType.h
#ifndef NODE_TYPE_H
#define NODE_TYPE_H

/**
 * @enum NodeType
 * @brief Specifies the type of node preprocessing when generating parameters.
 */
enum class NodeType {
    Equidistant,       /**< Equidistant node preprocessing */
    ChordalLength      /**< Cumulative chordal length node preprocessing */
};

#endif // NODE_TYPE_H