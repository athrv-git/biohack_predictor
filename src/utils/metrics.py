from sklearn.metrics import mean_absolute_error

def compute_mae(y_true, y_pred) -> float:
    """
    Compute Mean Absolute Error between true and predicted.
    """
    return mean_absolute_error(y_true, y_pred)
