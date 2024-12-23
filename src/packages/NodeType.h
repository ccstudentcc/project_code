// NodeType.h
#pragma once

/**
 * @enum NodeType
 * @brief Specifies the type of node preprocessing when generating parameters.
 */
enum class NodeType {
    Equidistant,       /**< Equidistant node preprocessing */
    ChordalLength      /**< Cumulative chordal length node preprocessing */
};

