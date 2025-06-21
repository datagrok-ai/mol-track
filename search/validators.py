"""
Additional validators for search functionality
"""

from typing import List, Dict, Any
import models


class SearchValidator:
    """Additional validation logic for search requests"""
    
    @staticmethod
    def validate_molecular_fields(request: models.SearchRequest) -> List[str]:
        """Validate molecular search specific requirements"""
        errors = []
        
        def check_conditions(filter_obj: models.LogicalNode):
            for condition in filter_obj.conditions:
                if isinstance(condition, models.AtomicCondition):
                    # Check for molecular operators
                    if condition.operator in ["IS SIMILAR", "CONTAINS", "IS CONTAINED"]:
                        # Must be used with structure field
                        if not condition.field.endswith('.structure'):
                            errors.append(f"Molecular operator '{condition.operator}' can only be used with structure fields")
                        
                        # IS SIMILAR requires threshold
                        if condition.operator == "IS SIMILAR" and condition.threshold is None:
                            errors.append(f"IS SIMILAR operator requires a threshold value")
                        
                        # Validate SMILES string
                        try:
                            from rdkit import Chem
                            mol = Chem.MolFromSmiles(condition.value)
                            if mol is None:
                                errors.append(f"Invalid SMILES string: {condition.value}")
                        except Exception:
                            errors.append(f"Could not validate SMILES string: {condition.value}")
                
                elif isinstance(condition, models.LogicalNode):
                    check_conditions(condition)
        
        if request.filter:
            check_conditions(request.filter)
        
        return errors
    
    @staticmethod
    def validate_cross_level_performance(request: models.SearchRequest) -> List[str]:
        """Check for potential performance issues with cross-level queries"""
        warnings = []
        
        if not request.filter:
            return warnings
        
        # Count distinct levels in conditions
        levels_used = set()
        
        def collect_levels(filter_obj: models.LogicalNode):
            for condition in filter_obj.conditions:
                if isinstance(condition, models.AtomicCondition):
                    level = condition.field.split('.')[0]
                    levels_used.add(level)
                elif isinstance(condition, models.LogicalNode):
                    collect_levels(condition)
        
        collect_levels(request.filter)
        
        # Warn about complex cross-level queries
        if len(levels_used) > 1:
            warnings.append("Cross-level query detected. Consider adding appropriate database indexes for optimal performance.")
        
        if len(levels_used) > 2:
            warnings.append("Complex multi-level query. This may have significant performance impact on large datasets.")
        
        return warnings
    
    @staticmethod
    def validate_output_consistency(request: models.SearchRequest) -> List[str]:
        """Validate that output fields are consistent with search level"""
        warnings = []
        
        # Check if output includes fields from other levels
        output_levels = set()
        for field in request.output:
            level = field.split('.')[0]
            output_levels.add(level)
        
        if len(output_levels) > 1:
            warnings.append("Output includes fields from multiple levels. Ensure proper JOIN relationships exist.")
        
        return warnings 