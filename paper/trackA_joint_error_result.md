# Track A joint label error simulation result

**Simulation. The candidate failed the joint robustness rule.** With known error structure and rates, joint correction passed every check. Under correlated error it brought null concordance from 0.274 for independent correction to 0.496. With underestimated rates, independent noise still moved null concordance to 0.547, outside the 0.03 tolerance. This separates the effects of knowing the correlation structure and knowing the marginal rates. [Saved results](../output/trackA_joint_error/result.json) | [Frozen plan](trackA_joint_error_plan.md)

The oracle had the simulator's error groups and rates. The next approach should use independent calibration or compare specimens after matching observation quality. Real data need independent region labels as well as coordinates.
