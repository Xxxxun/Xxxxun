## "Weight & Bias 网站使用心得"

简述wandb使用体验

---

#### Artifacts
wandb中每个run的artifacts保存后十分难以管理，强烈建议不要保存多个或不要在wandb上保存artifacts，否则难以批量删除。

#### wandb（心跳机制）
需要通知wandb关闭
```powershell
wandb.login(key="")
    try:
        wandb.init(
            project="PeSTo-training",
            config={
                "data_config": config_data,
                "model_config": config_model,
                "runtime_config": config_runtime
            },
            name=f"training_run_{config_runtime.get('experiment_name', 'default')}"
        )

        # train model
        train(config_data, config_model, config_runtime, './ckpt_pt')
    except KeyboardInterrupt:
        print("Interruption detected. Saving state...")
    except Exception as e:
        print(f"An error occurred: {e}")
        raise e
    finally:
        wandb.finish() 
        print("WandB run finished.")
```
```powershell
# 在终端运行，清理所有未完成的同步
wandb sync --clean
```
