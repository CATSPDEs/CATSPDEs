#pragma once
#include "AbstractQuadratureRule.hpp"

class MasterTriangleQuadratureRule 
	: public AbstractQuadratureRule<2, 3> {
	// singleton
	MasterTriangleQuadratureRule() : AbstractQuadratureRule {
		// (1) element (master triangle)
		{ { {0., 0.}, {1., 0.}, {0., 1.} } },
		// (2) nodes (n := polynomial order)
		{
			{ // n = 0
				{ 1/3., 1/3. }
			},
			{ // n = 1
				{ 1/3., 1/3. } 
			}, 
			{ // n = 2
				{ 1/6., 1/6. }, 
				{ 2/3., 1/6. }, 
				{ 1/6., 2/3. } 
			}, 
			{ // n = 3
				{ 1/3., 1/3. }, 
				{ 1/5., 1/5. }, 
				{ 1/5., 3/5. }, 
				{ 3/5., 1/5.} 
			}, 
			{ // n = 4
				{ .44594849091597, .44594849091597 },
				{ .44594849091597, .10810301816807 },
				{ .10810301816807, .44594849091597 },
				{ .09157621350977, .09157621350977 },
				{ .09157621350977, .81684757298046 },
				{ .81684757298046, .09157621350977 }
			}
		},
		{ // weights
			// n = 0
			{ 1/2. },
			// n = 1
			{ 1/2. },
			// n = 2
			{ 1/6., 1/6., 1/6. },
			// n = 3
			{ -27/96., 25/96., 25/96, 25/96 },
			// n = 4
			{ .111690794839005, .111690794839005, .111690794839005, .05497587182766, .05497587182766, .05497587182766 }
		}
	}
	{}
	MasterTriangleQuadratureRule(MasterTriangleQuadratureRule const &);
	MasterTriangleQuadratureRule& operator=(MasterTriangleQuadratureRule const &);
public:
	static auto& instance() {
		static MasterTriangleQuadratureRule single;
		return single;
	}
};