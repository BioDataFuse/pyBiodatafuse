"""Codes taken from DreamWalk repository: https://github.com/eugenebang/DREAMwalk."""

import logging
import pickle
from typing import List

import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, roc_auc_score
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier

from pyBiodatafuse.algorithms.DREAMwalk.utils import set_seed

logger = logging.getLogger(__name__)


def split_dataset(pairf: str, embeddingf: str, validr: float, testr: float, seed: int) -> tuple:
    """Split the dataset into train, validation, and test sets.

    :param pairf: Path to the pair file
    :param embeddingf: Path to the embedding file
    :param validr: Validation ratio
    :param testr: Test ratio
    :param seed: Random seed
    :return: Tuple of train, validation, and test sets
    """
    with open(embeddingf, "rb") as fin:
        embedding_dict = pickle.load(fin)

    xs, ys = [], []
    pair_df = pd.read_csv(pairf, sep="\t")

    skipped = set()
    for drug, dis, label in pair_df.values:
        try:
            drug_emb = embedding_dict[drug]
        except KeyError:
            skipped.add(drug)
            continue
        try:
            dis_emb = embedding_dict[dis]
        except KeyError:
            skipped.add(dis)
            continue

        xs.append(drug_emb - dis_emb)
        ys.append(int(label))

    logger.info(f"Skipped {len(skipped)} pairs due to missing embeddings.")
    logger.info(f"Total pairs: {len(pair_df)} | Remaining pairs: {len(pair_df) - len(xs)}")

    # dataset split
    x, y = {}, {}
    x["train"], x["test"], y["train"], y["test"] = train_test_split(
        xs, ys, test_size=testr, random_state=seed, stratify=ys
    )
    if validr > 0:
        x["train"], x["valid"], y["train"], y["valid"] = train_test_split(
            x["train"],
            y["train"],
            test_size=validr / (1 - testr),
            random_state=seed,
            stratify=y["train"],
        )
    else:
        x["valid"], y["valid"] = [], []

    return x, y


def return_scores(target_list: list, pred_list: List[float]) -> list:
    """Return metrics for the given target and prediction lists.

    :param target_list: list of target values
    :param pred_list: list of predicted values
    :return: list of metrics
    """
    metric_list = [accuracy_score, roc_auc_score, average_precision_score, f1_score]

    scores = []
    for metric in metric_list:
        scores.append(metric(target_list, pred_list))
    return scores


def predict_dda(
    embeddingf: str,
    pairf: str,
    modelf: str = "clf.pkl",
    seed: int = 42,
    validr: float = 0.1,
    testr: float = 0.1,
):
    """Predict drug-disease association.

    :param embeddingf: Path to the embedding file
    :param pairf: Path to the pair file
    :param modelf: Path to the model file
    :param seed: Random seed
    :param validr: Validation ratio
    :param testr: Test ratio
    """
    set_seed(seed)
    x, y = split_dataset(pairf, embeddingf, validr, testr, seed)

    clf = XGBClassifier(
        base_score=0.5,
        booster="gbtree",
        eval_metric="error",
        objective="binary:logistic",
        gamma=0,
        learning_rate=0.1,
        max_depth=6,
        n_estimators=500,
        tree_method="auto",
        min_child_weight=4,
        subsample=0.8,
        colsample_bytree=0.9,
        scale_pos_weight=1,
        max_delta_step=1,
        seed=seed,
    )

    clf.fit(x["train"], y["train"])

    preds = {}
    scores = {}
    for split in ["train", "valid", "test"]:
        preds[split] = clf.predict_proba(np.array(x[split]))[:, 1]
        scores[split] = return_scores(y[split], preds[split])
        logger.info(
            f"{split.upper():5} set | Acc: {scores[split][0]*100:.2f}% | AUROC: {scores[split][1]:.4f} | AUPR: {scores[split][2]:.4f} | F1-score: {scores[split][3]:.4f}"
        )

    with open(modelf, "wb") as fw:
        pickle.dump(clf, fw)
