## 🚀 What’s this PR do?

The PR must explain which problem it solves or which feature it introduces. This can also be done in an issue; if this is the case, then it's ok for PR to just have closes #whatever

---

### 📎 Related issues (optional)

Closes #123

## 🧠 Context

It would be helpful to point out if you solved the problem in an interesting way. For example, I've reordered the fields in a struct to achieve better struct packing, which reduced memory usage and improved CPU caching.
Especially if you reached for some unorthodox solution or used one that at first glance seems weird/wrong, it's helpful to point this out in the PR.  Example: Usually, I would use locking primitives here, but because of {this specific reason}, we can get away w/o locking.
If anything you think is worth pointing out or that the reviewer might find interesting or helpful, feel free to put it into the description. Example: I've decided to use RWLock instead of Mutex because Reader to Writer ratio is 5:1.

## ✅ Checklist

- [ ] I’ve tested this locally
- [ ] I’ve added relevant docs or comments
- [ ] I’ve updated or created tests if needed
- [ ] The branch with the feature is named `feature/<name>`